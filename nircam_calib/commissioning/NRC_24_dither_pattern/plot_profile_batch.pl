#!/usr/bin/perl -w
use strict;
use warnings;

my $car;
my $nplot = 0;
my $plots = 0;
my @lw_filter;
my @sw_filter;
my @simulations;
#
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
$nplot = run_batch(\@simulations, \@sw_filter, \@sw_sca, \@type);
print "SW plots: $nplot\n";
$plots = $plots+ $nplot;
#
$nplot = run_batch(\@simulations, \@lw_filter, \@lw_sca, \@type);
print "LW plots: $nplot\n";
$plots = $plots+ $nplot;
print "total plots: $plots\n";

exit(0);
#----------------------------------------------------------------------
sub run_batch{
    my ($sim_ref, $filter_ref, $sca_ref, $type_ref) = @_ ;
    my @simulations = @$sim_ref;
    my @filter      = @$filter_ref;
    my @sca         = @$sca_ref;
    my @type        = @$type_ref;
    my $command ;
    my $ii;
    my $jj;
    my $kk;
    my $ll;
    my $ntype;
    my $nplot = 0;
    
    for($ll=0 ; $ll <= $#sca; $ll++){
	for($kk =0; $kk <= $#simulations; $kk++) {
	    $ntype = $#type;
	    if($simulations[$kk] eq 'mirage'){
		$ntype = $ntype - 1;
	    }
	    for($jj = 0 ; $jj <= $#filter ; $jj++){
		for($ii = 0 ; $ii <= $ntype ; $ii++){
		    $command = join(' ','plot_profiles.py','--sca',$sca[$ll],'--sim',$simulations[$kk],'--filter',$filter[$jj],'--type',$type[$ii]);
		    print "$command\n";
		    system($command);
		    $nplot++;
		}
	    }
	}
    }
    return $nplot;
}
				    
				 
			    
			   
