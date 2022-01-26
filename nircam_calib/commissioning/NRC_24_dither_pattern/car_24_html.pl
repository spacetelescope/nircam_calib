#!/usr/bin/perl -w
use strict;
use warnings;

my $date ;
my $file;
my $filter;
my $header;
my $html;
my $key;
my $line;
my $nlines;
my $observation;
my $plot_list;
my $sca;
my $simulation;
my $table_header;
my $type;
#
my @data = ();
my %plot;
my %gslp;
my %gcal;
my %gi2d;
my %mcal;
my %mi2d;
#
$html = '/home/cnaw/commissioning/car_24_apt_01073/plots/table.html';
$plot_list = 'car_24_html.cat'; 
#
my $ndata = 0;
open(CAT,"<$plot_list") || die "cannot find $plot_list";
while(<CAT>) {
    chop $_;	
    ($sca, $simulation, $filter, $type, $observation, $file) =
	split(' ',$_);
    push(@data, $_);
    $ndata++;
}
close (CAT);

#
open(HTML,">$html") || die "cannot open $html";
print HTML "<html>\n";
print HTML "<head>\n";
print HTML "<style>\n";
print HTML "table,th,td{\n border: 3px solid black\nborder-collapse: collapse;}\n";
 print HTML "tr:nth-child(even) {background-color: #dddddd;}\n";
#print HTML "table,th,td{\n padding: 5px;}\n";

print HTML "th {text-align: center;}\n";
print HTML "td {text-align: center;}\n";

print HTML "</style>\n";
print HTML "</head>\n";
print HTML "<title>PSF Profiles CAR 24 APT 1073 </title>\n";
print HTML "<body>\n";
#print HTML "<center>\n";
#print HTML "<h1>Absolute Magnitude of the Sun for several filters</h1>\n";
#
print HTML '<table style="width:100%" class="table">',"\n";
#print HTML "<colgroup>\n";
#print HTML '<col class="grey" />',"\n";
#print HTML '<col class="red" />',"\n";
#print HTML '<col class="blue" />',"\n";
#print HTML "</colgroup>\n";
#
print HTML "<thead>\n";
$table_header =' <tr> <th>SCA</th> <th>Filter</th> <th>Simulation</th> <th>Type</th> <th>Observation</th><th>plot</th></tr>';
$table_header =' <tr> <th>Mirage cal</th> <th>Mirage i2d</th> <th>Guitarra cal</th> <th>Guitarra i2d</th> <th>Guitarra slp</th></tr>';
print HTML $table_header,"\n";
print HTML "</thead>\n";
#
print HTML "<tbody>\n";

for(my $ii = 0 ; $ii <= $#data ; $ii++) {
    ($sca, $simulation, $filter, $type, $observation, $file) =
	split(' ', $data[$ii]);
    $key = join('_',$observation,$sca);
    my @junk = split('\/',$file);
    $file = $junk[$#junk];
    if($simulation eq 'guitarra'){
	if($type eq 'slp') {$gslp{$key} = $file;}
	if($type eq 'cal') {$gcal{$key} = $file;}
	if($type eq 'i2d') {$gi2d{$key} = $file;}
    }
    
    if($simulation eq 'mirage'){
	if($type eq 'cal') {$mcal{$key} = $file;}
	if($type eq 'i2d') {$mi2d{$key} = $file;}
    }
}
my $style = 'style="width: 55vw; min-width: 110px;"';
$nlines = 0;
my $null = 0;
foreach $key (sort(keys(%mcal))){
    if($mcal{$key} eq 'None') {
	$null++;
	next;
    }
    $line = '<tr><td>';
    $line = join(' ',$line,'<img src=',$mcal{$key},'></td><td>');
    $line = join(' ',$line,'<img src=',$mi2d{$key},'></td><td>');
    $line = join(' ',$line,'<img src=',$gcal{$key},'></td><td>');
    $line = join(' ',$line,'<img src=',$gi2d{$key},'></td><td>');
    $line = join(' ',$line,'<img src=',$gslp{$key},'></td>');
    $line = join(' ',$line,'</tr>');
    print HTML $line,"\n";
    $nlines++;
}

print HTML "</table>\n";
print HTML "<hr>\n";
#
# Close HTML file
#
my @tm = localtime();
$date = sprintf(" %04d-%02d-%02d %02d:%02d:%02d MST", $tm[5]+1900, $tm[4]+1, $tm[3], $tm[2], $tm[1],$tm[0]);
$line = "<address>cnaw (at) as (dot) arizona (dot) edu</address>\n";
print HTML $line;
$line = "<!-- hhmts start -->\n";
print HTML $line;
$line = join(' ','Last modified:', $date);
print HTML $line,"\n";
$line = "<!-- hhmts end -->\n";
print HTML $line;
$line = "</body> </html>\n";
print HTML $line;
close(HTML);
print "nlines $nlines  ndata $ndata null $null\n";
