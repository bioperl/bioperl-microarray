use strict;
use vars qw($OK);

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't', '.';
    }

    eval {
           require Class::MakeMethods;
           require Class::MakeMethods::Emulator::MethodMaker;
           require enum;
    };

    if( $@ ){
      $OK = 0;
      use Test;
      plan tests => 1;
    } else {
      $OK = 1;
      use Test;    
      plan tests => 9;
    }
}

if(!$OK){ skip(1,1); exit }

use Bio::Expression::Microarray::Affymetrix::Array;
ok 1;
use Bio::Expression::Microarray::Affymetrix::Data;
ok 2;
$Bio::Expression::Microarray::Affymetrix::Array::DEBUG = 1;
ok 3;
$Bio::Expression::Microarray::Affymetrix::Data::DEBUG = 1;
ok 4;

use Bio::Expression::MicroarrayIO;
ok 5;

open(STREAM, 'cat ./eg/133a.cel ./eg/133a.cel |');

my $affx = Bio::Expression::MicroarrayIO->new(
#						-file     => './eg/95av2_combo.cel',
						-fh       => \*STREAM,
#						-template => './eg/95av2.cdf',
						-template => './eg/133a.cdf',
						-format   => 'affymetrix',
					   );
ok 6;
#$affx->verbose(1);
ok 7;
my $array = $affx->next_array;
ok 8;

my $array = $affx->next_array;
ok 9;
