use strict;

use Bio::Expression::Microarray::Affymetrix::Array;
use Bio::Expression::Microarray::Affymetrix::Data;
$Bio::Expression::Microarray::Affymetrix::Array::DEBUG = 1;
$Bio::Expression::Microarray::Affymetrix::Data::DEBUG = 1;
use Bio::Expression::MicroarrayIO;

my $affx = Bio::Expression::MicroarrayIO->new(
						-file     => './eg/133b.cel',
						-template => './eg/133b.cdf',
						-format   => 'affymetrix',
					   );
my $array = $affx->next_array;

open(T,'>133b.featuregroups');
foreach my $f ($array->each_featuregroup){
  print T $f->id, "\n";
}
close(T);

open(T,'>133b.qcfeaturegroups');
foreach my $f ($array->each_qcfeaturegroup){
  print T $f->id, "\n";
}
close(T);

