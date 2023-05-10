#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;

use Getopt::Euclid qw(:vars);
use Smart::Comments '###';

use Const::Fast;
use File::Basename;
use List::AllUtils qw(sum max mesh bundle_by partition_by);
use POSIX qw(log1p);
use Storable qw(lock_store lock_retrieve);

# should be fine considering small enough $N
const my $EPSILON  => 1e-12;


### Loading database from cache file: $ARGV_database
my $database = lock_retrieve $ARGV_database;
my ($theta, $binarize, $seq_n, $space, $utags, $M, $N, $S, $P, $details)
    = @{$database}{ qw(theta binarize seq_n space utags M N S P details) };

### database properties
### $seq_n
### $space
### $binarize
### $theta
### $M
### $N
### $S

### configuration
### $ARGV_split_tags
### $ARGV_prior_type
### $ARGV_max_read_n
### $ARGV_max_vote_n

# (optionally) handle multiple ontologies
my %tags_for                    # by ontology (or for all)
    = $ARGV_split_tags
    ? partition_by { (split '::')[0] } @$utags
    :                (     'all'     => $utags )
;
#### %tags_for
my @onts = keys %tags_for;
### @onts
my $ont_n = @onts;
### $ont_n


### Computing tag priors (log Pe)

my %num_for = map {                 # by tag
       $_ => $ARGV_prior_type eq 'uniform'
           ? 1                          # 1 tag  (equal weight)
           : $S->{$_}                   # N obs for one tag
} @$utags;
#### %num_for

my %sums_for = bundle_by {          # by ontology (or for all)
    $_[0] => $ARGV_prior_type eq 'uniform'
           ? scalar @{$_[1]}            # N tags (equal weight)
           : sum @{$S}{ @{$_[1]} }      # N obs for tags
} 2, %tags_for;
#### %sums_for

my %den_for  = bundle_by {          # by tag
    map { $_ => $sums_for{$_[0]} } @{$_[1]}
} 2, %tags_for;
#### %den_for

my %Pe = map { $_ => $num_for{$_} / $den_for{$_} } @$utags;     # by tag
### assert: abs($ont_n - sum values %Pe) < $EPSILON
   %Pe = map { $_ =>  log $Pe{$_}                } @$utags;     # log scale
### %Pe


# build indices
# because Pij are ArrayRef[ArrayRef]
my $i = 0;
my %ridx_for = map { $_->[2] => $i++ } @$details;
### assert: $i == $M
my $j = 0;
my %tidx_for = map { $_      => $j++ } @$utags;
### assert: $j == $N

#### %ridx_for
#### %tidx_for


my %posts_for;

for my $infile (@ARGV_infiles) {
    ### Processing read mapping file: $infile

    my %matches_for;
    my $read_n = 0;
    my $vote_n = 0;
    my $curr_read = q{};

    open my $in, '<', $infile;

    LINE:
    while (my $line = <$in>) {
        # skip header
        next LINE if $line =~ m/^@/xms;

        chomp $line;
        my ($read, undef, $repr, $sam_fields) = split "\t", $line, 4;

        # track unique read count...
        # ... and stop when enough reads
        if ($read ne $curr_read) {
            last LINE if defined $ARGV_max_read_n
                   && $read_n >= $ARGV_max_read_n;
            $read_n++;
            $curr_read = $read;
         }

        # skip unexact matches (for now)
        next LINE unless $sam_fields =~ m/\b MD:Z: \d+ \b/xms;

        # fetch representative seq for aligned read
        $vote_n++;
        push @{ $matches_for{$read} }, $repr;

        # stop when enough votes
        last LINE if defined $ARGV_max_vote_n
               && $vote_n >= $ARGV_max_vote_n;
    }

    #### %matches_for

    my %product_for = %Pe;          # by tag
    my $alig_n = 0;

    while (my ($read, $matches) = each %matches_for) {
        for my $repr_id (@$matches) {    # in case of multiple exact matches
            my $ridx = $ridx_for{$repr_id};
            $product_for{$_} += $P->[ $ridx ][ $tidx_for{$_} ] for @$utags;
        }   # sum of log P (joint P for sample seqs)
        $alig_n++;
    }
    #### %product_for

    ### $read_n
    ### $alig_n
    ### $vote_n

    # compute marginal likelihood (in log scale)
    # (optionally) handle multiple ontologies
    my %sums_for = bundle_by {      # by ontology
        $_[0] => lse_all( @product_for{ @{$_[1]} } )
    } 2, %tags_for;
    #### %sums_for
    my %sum_for  = bundle_by {      # by tag
        map { $_ => $sums_for{$_[0]} } @{$_[1]}
    } 2, %tags_for;
    #### %sum_for

    # compute log posteriors using Bayes' theorem
    my %post_for = bundle_by {      # by tag
        $_[0] => $_[1] - $sum_for{$_[0]}
    } 2, %product_for;
    #### %post_for

    # store log posteriors for sample
    my ($sample) = fileparse( $infile, qr/\.sam$/xms );
    $posts_for{$sample} = {
        read_n => $read_n,
        alig_n => $alig_n,
        vote_n => $vote_n,
        map { $_ => $post_for{$_} // 'NA' } @$utags     # no NA if priors
    };

}


### Storing posterior probs to TSV outfile: $ARGV_outfile

my @samples = sort keys %posts_for;
my @columns = ( qw(read_n alig_n vote_n), @$utags );

open my $out, '>', $ARGV_outfile;
say {$out} join "\t", ('sample', @columns);

for my $sample (@samples) {
    say {$out} join "\t", $sample,
        map { $posts_for{$sample}{$_} } @columns;
}


# from: https://stats.stackexchange.com/questions/379335/adding-very-small-probabilities-how-to-compute
sub lse {
    my $l1 = shift;
    my $l2 = shift;
    return max($l1, $l2) + log1p(exp(-abs($l1-$l2)));
}

sub lse_all {
    my $l1 = shift;
    while (my $l2 = shift) {
        $l1 = lse($l1, $l2);
    }
    return $l1;
}

__END__

=head1 USAGE

    comp_tag_probs_in_samples.pl <infiles> --database=<file> --outfile=<file> \
         [optional arguments]

=head1 REQUIRED ARGUMENTS

=over

=item <infiles>

Path to input SAM files [repeatable argument].

=for Euclid:
    infiles.type: readable
    repeatable

=item --database=<file>

Path to input database (cache) file.

=for Euclid:
    file.type: readable

=item --outfile=<file>

Path to TSV outfile with posterior probs for tags given sample seqs.

=for Euclid:
    file.type: writable

=back

=head1 OPTIONAL ARGUMENTS

=over

=item --prior-type=<str>

Type of tag priors to use [default: str.default]. The following types are
available:

    - database (priors computed from tag distribution in database)
    - uniform  (flat priors computed as 1/N)

=for Euclid:
    str.type:    /database|uniform/
    str.type.error: <str> must be one of database or uniform (not str)
    str.default: 'uniform'

=item --split-tags

Consider tags as belonging to distinct ontologies [default: no].
Ontology-aware tags are encoded as C<Ont::Tag>. Enabling this option also
affects the way tag priors are computed.

=item --max-read-n=<n>

Maximum number of input reads (aligned or not) to process [default: n.default].

=for Euclid:
    n.type: +number

=item --max-vote-n=<n>

Maximum number of votes (= read mappings) to process [default: n.default].
This number can grow faster than the number of aligned reads when multiple
mappings per read are allowed in the SAM file(s).

=for Euclid:
    n.type: +number

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back
