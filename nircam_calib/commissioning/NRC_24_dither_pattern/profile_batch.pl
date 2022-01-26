#!/usr/bin/perl -w
use strict;
use warnings;

my $batch_file = 'profile.batch';
my $car_root;
my $command;
my $cutoff = 0.0;
my $file ;
my $host;
my $overwrite = 1;
my $profile;
my $target_dir ="." ;
my $type;
my $debug        = 'false';
my $show_plot    = 'false' ;
my $use_sep      = 'true' ;
my @list ;
if($#ARGV <= -1) {
    print "./run_profile_batch.pl mirage cal\n";
    print "./run_profile_batch.pl guitarra cal parallel\n";
    print "./run_profile_batch.pl templates parallel\n";
    exit(0);
}
my $sim = $ARGV[0];
if($sim ne 'templates') {
    $type = $ARGV[1];
}
my $parallel = 0;
for(my $ii = 2; $ii <= $#ARGV; $ii++){
    if(lc($ARGV[$ii]) eq 'parallel') {
	$parallel = 1;
    }
}
#
$host = $ENV{HOST};
print "host is $host\n";
$car_root = $ENV{CAR_ROOT};
if($host eq 'ema.as.arizona.edu'){
    $car_root = '/home/cnaw/commissioning/car_24_apt_01073/';
    if($sim eq 'mirage') {
	$target_dir = $car_root.'mirage/analysis/';
    }
    if($sim eq 'guitarra') {
	$target_dir = $car_root.'guitarra/analysis/';
    }
    if($sim eq 'templates') {
	$target_dir = '/home/cnaw/python/commissioning/templates/';
    }
}
if($host eq 'orange.as.arizona.edu'){
    $car_root = '/data1/car_24_apt_01073/';
    if($sim eq 'mirage') {
	$target_dir = $car_root.'mirage/analysis/';
    }
    if($sim eq 'guitarra') {
	$target_dir = $car_root.'/guitarra/analysis/';
    }
    if($sim eq 'templates') {
	$target_dir = '/home/cnaw/python/commissioning/templates/';
    }
}
print "car_root is defined as $car_root\n";

#@list = `ls $target_dir*cal*psf.fits | grep F070W`;
#@list = `ls $target_dir*cal*psf.fits`;
#@list = `ls $target_dir*slp*psf.fits`;
#@list = `ls $target_dir*cal*psf.fits`;
if($sim eq 'templates') {
    my $string = $target_dir.'PSF_NIRCam_F*W*.fits';
    @list = `ls $string | grep -v prof |grep -v mask`;
    $use_sep = 'false';
    $use_sep = 'true';
} else{
    @list = `ls $target_dir*$type*psf.fits`;
}
if($parallel == 1) {
    open(BATCH,">$batch_file") || die "cannot open $batch_file";
}
for(my $ii = 0 ; $ii <= $#list ; $ii++) {
    $file = $list[$ii];
    $file =~ s/\n//g;
    $profile = $file;
    $profile =~ s/_psf.fits/_psf_prof.fits/;
    if(-e $profile && $overwrite == 0) {
	print "$profile exists; skipping\n";
	next ;
    } else {
	$command = join(' ','single_psf.py','--file', $file);
	if(lc($use_sep eq 'true')){
	    $command = join(' ',$command, '--use_sep');
	}
	if(lc($show_plot eq 'true')){
	    $command = join(' ',$command, '--plot');
	}
	if(lc($debug eq 'true')){
	    $command = join(' ',$command, '--debug');
	}
	if($parallel == 1) {
	    print BATCH $command,"\n";
	} else {
	    print "$command\n";
	    system($command);
	}
    }
#    exit(0);
}
if($parallel == 1){
    my $line = join(' ','./parallel.pl ', $batch_file);
    print "$line\n";
}
    
   
