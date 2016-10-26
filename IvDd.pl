#!/usr/local/bin/perl

#####################################################
# Estimating putative insertinal variations of HML-2 from DGV archive.
# this.pl -d DGV_file -t LTR_list.bed -i ORF_list.bed [other_options]
# Results are printed on the standard output.
# Below tools and data are required in this script.
# [1] bedtools (intersectBed)	# the latest version recommended. -F option is utilized for intersectBed.
# [2] DGV dataset (.txt)	# from DGV Variants (20 columns)
# [3] LTR list (.bed)		# from RepeatMasker (6 columns)
# [4] int list (.bed)		# from RepeatMasker (6 columns)
# Table captions of the result: <6 captions of [3]> <Chr, Start, End, ID, CNV, Type and Study of [2]> <6 captions of [4]>
#####################################################

use strict;
use warnings;
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);
use Pod::Usage;

my $path_intersectBed;	# path of intersectBed
my $db_origin;		# DGV dataset (.txt)
my $target;		# LTR list (.bed)
my $int;		# Internal proviral sequence list (.bed) 
my $delen = 1500;	# Maximum of deletion length parameter
my $ovlp = 1;		# intersectBed -f option
my $gap_LTR_int_max = 100;
my %opt = ();
GetOptions(\%opt, 'b=s' => \$path_intersectBed,'f=s' => \$ovlp, 'd=s' => \$db_origin, 't=s' => \$target, 'del=i' => \$delen, 'i=s' => \$int, 'help|h');
pod2usage(1) if ($opt{help} or $opt{h});

die "ERROR: -b, -d, -t and -i options required!\n" if(!defined $path_intersectBed or !defined $db_origin or !defined $target or !defined $int);
my $db = 'db.tmp';
system("cat $db_origin | perl -lane '\$chr=\"chr\".\$F[1];print join \"\t\", (\$chr,\@F[2,3],\$F[0],\@F[4..6])' | grep -e deletion -e loss >$db");

# solo_LTR
my $tmp = 'tmp.bed';
system("$path_intersectBed -a $target -b $db -f $ovlp -wao | grep -v '\\.'>$tmp");
my $int_col_num = `head -1 $int | perl -lane '\$count=1;while(\$_=~/\t/g){\$count++};print \$count'`;
my @empty = split //, ('*' x $int_col_num);
my $data;
open my $t, '<', $tmp or die $!;
while(my $line = <$t>) {
	chomp $line;
	my @content = split /\t/, $line;
	next if($content[7] < 0);
	
	my ($ltr_s, $ltr_e, $del_s, $del_e) = @content[1,2,7,8];
	if($ltr_s < $del_s or $ltr_e > $del_e) {
		next;
	}
	
	my $del_len = $del_e - $del_s + 1;
	next if($del_len >= $delen);

	my $key = join ":", @content[0..2];
	if(exists $data->{$key}) {
		$del_len = (split /\t/, $data->{$key})[-1] if((split /\t/, $data->{$key})[-1] < $del_len);
	}
	$data->{$key} = join "\t", (@content[0..12], @empty, $del_len);
}
close $t;

# provirus
my $tmp2 = 'tmp2.bed';
my $tmp3 = 'tmp3.bed';
system("$path_intersectBed -a $db -b $int -F 1 -wao | grep -v '\\.' >$tmp2");
system("$path_intersectBed -a $target -b $tmp2 -f $ovlp -wao | grep -v '\\.' >$tmp3");

open my $tt, '<', $tmp3 or die $!;
while(my $line = <$tt>) {
	chomp $line;
	my @content = split /\t/, $line;
	next if($content[7] < 0);	

	my ($ltr_s, $ltr_e, $del_s, $del_e, $int_s, $int_e) = @content[1,2,7,8,14,15];
	my $ltr_int_len = 0;
	if(($ltr_e <= $int_s and $int_s - $ltr_e > $gap_LTR_int_max) or ($ltr_s >= $int_e and $ltr_s - $int_e > $gap_LTR_int_max)) {
		next;
	} elsif($ltr_e <= $int_s) {
		$ltr_int_len = $int_e - $ltr_s + 1;
	} elsif($ltr_s >= $int_e) {
		$ltr_int_len = $ltr_e - $int_s + 1;
	}
	
	my $del_len = $del_e - $del_s + 1;
	my $corrected_del_len = $del_len - $ltr_int_len;
	next if($corrected_del_len >= $delen);
	
	my $key = join ":", @content[0..2];
	if(exists $data->{$key}) {
		my $check = (split /\t/, $data->{$key})[19];
		(@content[0..18], $corrected_del_len) = (split /\t/, $data->{$key}) if(defined $check and $check < $corrected_del_len);
	}
	$data->{$key} = join "\t", (@content[0..18], $corrected_del_len);
}

for my $v (sort {
	my $aa = (split /\t/, $a)[0]; $aa =~ s/chr//; $aa = 23 if($aa eq 'X'); $aa = 24 if($aa eq 'Y');
	my $bb = (split /\t/, $b)[0]; $bb =~ s/chr//; $bb = 23 if($bb eq 'X'); $bb = 24 if($bb eq 'Y');
	$aa <=> $bb
	} sort {
	my $aaa = (split /\t/, $a)[1];
	my $bbb = (split /\t/, $b)[1];
	$aaa <=> $bbb
	}values %{$data}) {
	my @splitv = split /\t/, $v;
	print join "\t", @splitv[0..18], "\n";
}

close $tt;

system("rm $db $tmp $tmp2 $tmp3");
exit;

__END__
 
=pod
 
=head1 SYNOPSIS
 
B<perl IvDd.pl -d dgv.txt -t RM_LTR_list.bed -i RM_int_list.bed [other_options]>
 
Options: [-h|--help|-b*|-d*|-f|-t*|-i*|--del; *required]

Output: Standard output

=head1 OPTIONS
 
=over 8
 
=item B<-h|--help>
 
print help.

=item B<-b>
 
/path/intersectBed

=item B<-d>
 
DGV data (.txt)

=item B<--del>
 
Maximum of deletion length parameter (default 1500)

=item B<-f>
 
intersectBed -f option (default 0.9)

=item B<-t|-i>
 
LTR_list.bed (-t) and int_list.bed (-i)

=back

=head1 DESCRIPTION

=cut
