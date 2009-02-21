# $Id$
# BioPerl module for Bio::MicroarrayIO::cel
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::Affymetrix::Data - Affymetrix CEL file.

=head1 SYNOPSIS

You should not be using this module directly.  Try
Bio::Expression::MicroarrayIO instead.

=head1 DESCRIPTION

Bio::Expression::Microarray::Affymetrix::Data just holds a reference to
a big matrix of values corresponding to the spots on an Affy array.  It
also knows how to poke data into the matrix slots.

=head1 FEEDBACK

Direct feedback to E<lt>allenday@ucla.eduE<gt> or to the Bioperl mailing list (see below).

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 
 
Please direct usage questions or support issues to the mailing list:
  
L<bioperl-l@bioperl.org>
  
rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::Expression::Microarray::Affymetrix::Data;

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Expression::Microarray::Affymetrix::Feature;
use IO::File;

use base qw(Bio::Root::Root Bio::Root::IO);
use vars qw($DEBUG);

use Class::MakeMethods::Emulator::MethodMaker
  get_set => [qw( array mode )],
  new     => 'new',
;

sub load_data {
  my($self,$line) = @_;

  my($key,$value) = (undef,undef);
  if(my($try) = $line =~ /^\[(.+)\]/){
	$self->mode($try);
	return;
  }

  ($key,$value) = $line =~ /^(.+?)=(.+)$/;

  if($self->mode eq 'CEL'){
	$self->array->cel($self->array->cel . $line) if $key;

	$self->array->version($value) if $key eq 'Version';
  } elsif($self->mode eq 'HEADER'){
	$self->array->header($self->array->header . $line) if $key;

	$self->array->algorithm($value) if $key eq 'Algorithm';
	$self->array->algorithm_parameters($value) if $key eq 'AlgorithmParameters';
	$self->array->dat_header($value) if $key eq 'DatHeader';
  }
  elsif($self->mode eq 'INTENSITY'){
	if($key){
	  $self->array->intensity($self->array->intensity . $line);
	  return;
	}
	
	$line =~ s/\s*(.+)\s*/$1/; #clean up spaces;
	my @row = split /\s+/, $line;
	
	print STDERR sprintf("%-60s\r",sprintf("%3d %3d is %4.1f",$row[0],$row[1],$row[2])) if $DEBUG;
	
	return unless defined $row[1] and defined $row[2];
	my $feature = $self->array->matrix($row[0],$row[1]);
	
	if(defined $feature){
	  $$feature->value($row[2]);
	  $$feature->standard_deviation($row[3]);
	  $$feature->sample_count($row[4]);
	} else {
	  #warn "$row[0] $row[1] INTENSITY" if $DEBUG;

	  #There are some spots on the array (CEL file) that do not have
	  #a corresponding entry in the CDF file.  It is assumed that these
	  #are quality control features of some sort.  We opt to create feature
	  #objects for them, because their values are necessary to recreate a
	  #CEL file.

	  $feature = Bio::Expression::Microarray::Affymetrix::Feature->new(
																	   x => $row[0],
																	   y => $row[1]
																	  );
	  $feature->is_singleton(1);
	  $feature->value($row[2]);
	  $feature->standard_deviation($row[3]);
	  $feature->sample_count($row[4]);

	  $self->array->matrix($row[0],$row[1],\$feature);
	}
  }
  elsif($self->mode eq 'MASKS'){
	$self->array->masks($self->array->masks . $line) and return if $key;

	my @row = split /\t/, $line;
	return if $row[0] == 0 and $row[1] == 0;		# WHY???
	my $feature = $self->array->matrix($row[0],$row[1]);
	$$feature->is_masked(1) if ref($feature);		# WHY???
  }
  elsif($self->mode eq 'OUTLIERS'){
	$self->array->outliers($self->array->outliers . $line) and return if $key;
	
	my @row = split /\t/, $line;
	my $feature = $self->array->matrix($row[0],$row[1]);	
	$$feature->is_outlier(1) if ref($feature);		# WHY???
  }
  elsif($self->mode eq 'MODIFIED'){
	$self->array->modified($self->array->modified . $line) and return if $key;
	
	my @row = split /\t/, $line;
	my $feature = $self->array->matrix($row[0],$row[1]);
	$$feature->is_modified(1) if ref($feature);		# WHY???
  }
}

1;
