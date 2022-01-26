#!/usr/bin/perl -w
use strict;
use warnings;

my $car;
my $file_ref;
my $html_cat;
my $key;
my $line;
my $matrix_ref;
my $nplot = 0;
my $number_obs = 7;
my $plots = 0;

my @matrix;
my @lw_filter;
my @sw_filter;
my @simulations;

my %table;
#
$html_cat = 'car_24_html.cat';
my @type    = ('cal','i2d','slp');
my @lw_sca  = ('a5','b5');
my @sw_sca  = ('a1','a2','a3','a4','b1','b2','b3','b4');

$car = 24;

if($car == 24){
    @sw_filter = ('F070W','F115W','F150W');
    @lw_filter = ('F277W');
    @simulations = ('guitarra','mirage');
}
#
($nplot, $matrix_ref, $file_ref) = run_batch(\@simulations, \@sw_filter, \@sw_sca, \@type, $number_obs);
print "SW plots: $nplot\n";
$plots = $plots+ $nplot;
#
($nplot, $matrix_ref, $file_ref) = run_batch(\@simulations, \@lw_filter, \@lw_sca, \@type, $number_obs);
print "LW plots: $nplot\n";
$plots = $plots+ $nplot;
#

@matrix = @$matrix_ref;
open(CAT,">$html_cat") || die "cannot open $html_cat";
for(my $ii =0 ; $ii <= $#matrix ; $ii++){
    my $file = $file_ref->[$ii];
    if($file =~ m/None/) { next;}
    $key  = $matrix[$ii];
#    $table{$key} = $file;
#    ($sca, $simulation, $filter, $type, $observation) =
#	split(' ',$key);
    $line = join(' ', $key, $file);
    print CAT $line,"\n";
}
print "total plots: $plots\n";

exit(0);
#----------------------------------------------------------------------
sub run_batch{
    my ($sim_ref, $filter_ref, $sca_ref, $type_ref, $number_obs) = @_ ;
    my @simulations = @$sim_ref;
    my @filter      = @$filter_ref;
    my @sca         = @$sca_ref;
    my @type        = @$type_ref;
    my $command ;
    my $hh;
    my $ii;
    my $jj;
    my $kk;
    my $ll;
    my $ntype;
    my $nplot = 0;
    my $observation;
    my @files = ();
    my @matrix = ();
    my $output = 'temp.list';
    unlink $output;
    $command = join(' ','touch',$output);
    system($command);
    #   sca -> simulation -> type -> filter -> observation
    #   or
    #   sca -> filter -> observation  -> simulation -> type ?
    for($ll=0 ; $ll <= $#sca; $ll++){
	for($kk =0; $kk <= $#simulations; $kk++) {
	    $ntype = $#type;
	    if($simulations[$kk] eq 'mirage'){
		$ntype = $ntype - 1;
	    }
	    for($jj = 0 ; $jj <= $#filter ; $jj++){
		for($ii = 0 ; $ii <= $ntype ; $ii++){
		    for($hh = 1 ; $hh <= $number_obs ; $hh++) {
			$observation = sprintf("%03d", $hh);
			$command = join(' ','plot_profiles.py','--sca',$sca[$ll],'--sim',$simulations[$kk],'--filter',$filter[$jj],'--type',$type[$ii],'--obs', $observation, ">>",$output);
			print "$command\n";
			system($command);
			$nplot++;
			my $line = join(' ',$sca[$ll],$simulations[$kk], $filter[$jj], $type[$ii],$observation);
			push(@matrix, $line);
		    }
		}
	    }
	}
    }
    open(OUT,"<$output") || die "cannot open $output";
    my $nn = -1;
    while(<OUT>) {
	chop $_;
	push(@files,$_);
	$nn++;
	print "$matrix[$nn], $_\n";
    }
    close(OUT);
    return $nplot, \@matrix, \@files;
}
