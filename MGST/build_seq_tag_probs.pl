#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;

use Getopt::Euclid qw(:vars);
use Smart::Comments '###';

use Const::Fast;
use List::AllUtils qw(all count_by);
use Storable qw(lock_store lock_retrieve);

use Bio::MUST::Core::Utils qw(:filenames);
use aliased 'Bio::FastParsers::CdHit';

const my $STORABLE => '.storable';
const my $EPSILON  => 1e-10;            # needed for large $M


### Processing TSV seq-tags file: $ARGV_seq_tags
my $tags_for;
my $tags_for_store = append_suffix($ARGV_seq_tags, $STORABLE);

if (-e $tags_for_store) {
    ### Loading from cache file: $tags_for_store
    $tags_for = lock_retrieve $tags_for_store;
}
else {
    open my $tag_in, '<', $ARGV_seq_tags;
    while (my $line = <$tag_in>) {      ### Processing seqs |===[%]
        chomp $line;
        my ($acc, $tag_str) = split "\t", $line;
        my @tags = split ',', $tag_str;
        $tags_for->{$acc} = \@tags;
    }
    ### Storing to cache file: $tags_for_store
    lock_store $tags_for, $tags_for_store;
}
##### $tags_for


### Processing CD-HIT cluster file: $ARGV_clusters
my $report;
my $report_store = append_suffix($ARGV_clusters, $STORABLE);

if (-e $report_store) {
    ### Loading from cache file: $report_store
    $report = lock_retrieve $report_store;
}
else {
    $report = CdHit->new( file => $ARGV_clusters );
    ### Storing to cache file: $report_store
    lock_store $report, $report_store;
}
my @representatives = $report->all_representatives;

### Collecting observations (bij) for representative seqs
my @details;
my @B;
my %S;

### $ARGV_binarize

my $clus_id = 0;
my $seq_n = 0;
my $space = 0;

for my $repr_id (@representatives) {    ### Collecting |===[%]
    #### $clus_id
    #### $repr_id
    my $members = $report->members_for($repr_id);
    my $clus_sz = @$members + 1;
    #### $members
    #### $clus_sz
    $seq_n += $clus_sz;
    my @tags = map { @{ $tags_for->{$_} } } $repr_id, @$members;
    #### @tags
    my $tag_n = @tags;
    $space += $tag_n;
    push @details, [ $clus_id++, $clus_sz, $repr_id, $tag_n ];
    my %count_for = count_by { $_ } @tags;
    if ($ARGV_binarize) {
        %count_for = map { $_ => 1 } keys %count_for;
    }
    #### %count_for
    push @B, \%count_for;
    $S{$_} += $count_for{$_} for keys %count_for;
    #### %S
}
#### @B

my $M = @B;
my @utags = sort keys %S;
my $N = @utags;

### $seq_n
### $space

### $M
### $N

### %S


### Computing and storing cond probs (log Pij) to TSV outfile: $ARGV_outfile
my @P;

my $theta = 1 - $ARGV_error_rate;
### $ARGV_error_rate
### $theta

open my $out, '>', $ARGV_outfile;
say {$out} join "\t", qw(clus_id clus_sz repr_id tag_n),                @utags;
say {$out} join "\t",    $M,     $seq_n, 'S',    $space, map { $S{$_} } @utags;

my @sums;

for (my $i = 0; $i < @B; $i++) {        ### Computing |===[%]
    my @Pi = map {
        ( (1 - $theta) / $M ) + ( $theta * ($B[$i]{$_} // 0) / $S{$_} )
    } @utags;
    $sums[$_] += $Pi[$_] for 0..$#Pi;
    @Pi = map { log $_ } @Pi;
    say {$out} join "\t", @{ $details[$i] }, @Pi;
    push @P, \@Pi;
}
#### @sums
### assert: @P == $M
### assert: @sums == $N
### assert: all { abs(1 - $_ < $EPSILON) } @sums

my $database = {
    theta    => $theta,
    binarize => $ARGV_binarize,
    seq_n    => $seq_n,
    space    => $space,
    utags    => \@utags,
    details  => \@details,
    M  => $M,
    N  => $N,
    S  => \%S,
    P  => \@P,
};

my $database_store
    = insert_suffix( change_suffix($ARGV_outfile, $STORABLE), '-database' );
### Storing to cache file: $database_store
lock_store $database, $database_store;


__END__

=head1 USAGE

    build_seq_tag_probs.pl --seq-tags=<file> --clusters=<file> \
        --outfile=<file> [optional arguments]

=head1 REQUIRED ARGUMENTS

=over

=item --seq-tags=<file>

Path to input file with sequence => tag(s) pairs in TSV format.

=for Euclid:
    file.type: readable

=item --clusters=<file>

Path to input file with cluster => sequences(s) pairs in CD-HIT format.

=for Euclid:
    file.type: readable

=item --outfile=<file>

Path to TSV outfile (and database) with conditional probs for sequences given
tags.

=for Euclid:
    file.type: writable

=back

=head1 OPTIONAL ARGUMENTS

=over

=item --error-rate=<n>

Error rate (0 to 1) to be applied to observations when computing conditional
probs [default: n.default]. Confidence (theta) will be derived as 1 -
error-rate.

=for Euclid:
    n.type:    +number
    n.default: 1e-9

=item --binarize

Build a purely binary database of conditional probs (still taking into account
--theta) [default: no]. Otherwise, tags attached to multiple sequences of a
cluster are given more weight.

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back
