# $Id$
# BioPerl module for Bio::Expression::MicroarrayIO::affymetrix
#
# Copyright Allen Day <allenday@ucla.edu>, Stanley Nelson
# <snelson@ucla.edu>.  Human Genetics, UCLA Medical School,
# University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::MicroarrayIO::affymetrix - Work with an Affymetrix
array.

=head1 SYNOPSIS

  $stream  = Bio::Expression::MicroarrayIO->new(
		'-file'     => "my.cel",
		'-template' => "my.cdf",
		'-format'   => "affymetrix",
						);

=head1 DESCRIPTION

Bio::Expression::MicroarrayIO::affymetrix parses Affymetrix CEL files in
the context of a CDF file.

=head1 FEEDBACK

Direct feedback to E<lt>allenday@ucla.eduE<gt> or to the Bioperl mailing list (see below).

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

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
package Bio::Expression::MicroarrayIO::affymetrix;

use strict;
use Bio::Root::Root;
use Bio::Expression::MicroarrayIO;
use Bio::Expression::Microarray::Affymetrix::CelArray;
use Bio::Expression::Microarray::Affymetrix::ArrayDesign;
use IO::File;

use base qw(Bio::Root::Root Bio::Expression::MicroarrayIO);
use vars qw($started);

use constant CRLF => "\015\012";

=head2 new

 Title   : new
 Usage   : Bio::Expression::MicroarrayIO::affymetrix->new(
											  -file     => 'path/to/filename',
											  -template => 'path/to/template')
 Comments: You should probably not be instantiating this module directly.
           Use MicroarrayIO instead.
 Args    : -file     => filename
           -format   => format
           -template => template

=cut

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;

    my $self = bless {}, $class;

    $self->_initialize(@args);
    return $self;
}

=head2 _initialize

 Title   : _initialize
 Usage   : $affy->_initialize(
							  -file     => 'path/to/filename',
							  -template => 'path/to/template'
							 );
 Function: Loads up a template module that will be used to
           by next_array();
 Returns : nothing.
 Args    : -file     => filename
           -format   => format
           -template => template

=cut

sub _initialize{
  my ($self,@args) = @_;
  $self->SUPER::_initialize(@args);
  $self->_initialize_io(@args);

  my %param = @args;
  @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
  $self->verbose($param{-verbose});

  if($param{-file} or $param{-fh}){
	print STDERR "loading array template and data...\n" if $self->verbose;

	$self->datafile($param{-file});
	$self->templatefile($param{-template});

	print STDERR "array template loaded!...\n" if $self->verbose;
  }
}

=head2 templatefile

 Title   : templatefile
 Usage   : $affy->templatefile('path/to/template');
           $affy->templatefile();
 Function: get/set the the location of a template file
 Returns : path to a template file
 Args    : optional path to a template file

=cut

sub templatefile {
  my($self,$val) = @_;
  $self->{templatefile} = $val if $val;
  return $self->{templatefile};
}

=head2 datafile

 Title   : datafile
 Usage   : $affy->datafile('path/to/datafile');
           $affy->datafile();
 Function: get/set the the location of a data file
 Returns : path to a data file
 Args    : optional path to a data file

=cut

sub datafile {
  my($self,$val) = @_;
  $self->{datafile} = $val if $val;
  return $self->{datafile};
}

=head2 array

 Title   : array
 Usage   : $affy->array($template);
           $affy->array();
 Comments: You probably should not be using this method
           to set the array template object.  Use load_array() instead.
 Function: get/set the the array template object
 Returns : a Bio::Expression::Microarray::Affymetrix::ArrayDesign object
 Args    : optional Bio::Expression::Microarray::Affymetrix::ArrayDesign object

=cut

sub array {
  my($self,$val) = @_;
  $self->{array} = $val if $val;
  return $self->{array};
}

=head2 load_array

 Title   : load_array
 Usage   : $affy->load_array($template);
           $affy->load_array();
 Function: cause a Bio::Expression::Microarray::Affymetrix::ArrayDesign
           object to be created using $affy->templatefile().
 Returns : a  Bio::Expression::Microarray::Affymetrix::ArrayDesign object
 Args    : optional path to a template file, which is stored in
           $affy->templatefile before the Template object is created.

=cut

sub load_array {
  my($self,$arg) = @_;

  return $self->array if defined($self->array);

  warn "loading array..." if $self->verbose;
  $self->templatefile($arg) if defined $arg;

  my $array = Bio::Expression::Microarray::Affymetrix::ArrayDesign->new(
				  -file => $self->templatefile,
									 );
  my $t;
  open($t,$self->templatefile) or $self->throw("Couldn't open templatefile: ".$self->templatefile.": $!");
  while(<$t>){
	$array->load_data($_);
  }
  close($t);

  $self->array($array);
  return $self->array;
}

=head2 next_array

 Title   : next_array
 Usage   : $affy->next_array();
 Function: reads an Affymetrix data record from $affy->datafile
 Returns : a loaded array object
 Args    : none

=cut

sub next_array {
  my $self = shift;
  $self->load_array();
  $self->array->destroy_features;
  print STDERR "loading data...\n" if $self->verbose;

  #create an object for parsing the array data;
  my $data = new Bio::Expression::Microarray::Affymetrix::CelArray;

  $self->array->clear_cel();
  $self->array->clear_header();
  $self->array->clear_modified();
  $self->array->clear_intensity();
  $self->array->clear_masks();
  $self->array->clear_outliers();

  while( defined( $_ = $self->_readline(-raw=>1) ) ){
	#if we see the beginning of a new file
	if($_ =~ m!\[CEL\]!){
	  print STDERR "new CEL\n" if $self->verbose;
	  $data->array($self->array);

	  #if we have seen another CEL before, _pushback() the beginning
	  #of this new CEL, for the next call to next_array() to handle,
	  #and return the preceding CEL.
	  if($started){
		$self->_pushback($_);
		undef $started;
		return $self->array;
	  }
	  $started = 1;
	}
	next unless $started;

	#slurp up the data, line by line
	$data->load_data($_);
  }

  print STDERR "loaded data!\n" if $self->verbose;

  #and return the array object, complete with loaded data!
  return $self->array;
}

=head2 write_array

 Title   : write_array
 Usage   : $affy->write_array($array);
 Function: write an Affymetrix data record using $array
 Returns : nothing.  prints a lot of text to the MicroarrayIO stream.
 Args    : A Bio::Expression::MicroarrayI compliant object.

=cut

sub write_array {
  my($self,$array) = @_;

  my @masks;
  my @modified;
  my @outliers;

  if( !defined $array ) {
	$self->throw("Attempting to write with no array!");
  }

  if( ! ref $array || ! $array->isa('Bio::Expression::MicroarrayI') ) {
	$self->warn(" $array is not MicroarrayI compliant. Dump may fail!");
  }

  $self->_print("[CEL]"       . CRLF . $array->cel       . CRLF);
  $self->_print("[HEADER]"    . CRLF . $array->header    . CRLF);
  $self->_print("[INTENSITY]" . CRLF . $array->intensity);

  my $i = -1;
  foreach my $x ( @{ $array->matrix } ){
	$i++;
	next unless $x;
	my $j = 0-1;
	foreach my $y ( @$x ){
	  $j++;
	  next unless $y;

	  #$$y is a probe object
	  $self->_print(join "\t", (sprintf("%3d",$j),
						sprintf("%3d",$i),
						sprintf("%.1f",$$y->value),
						sprintf("%.1f",$$y->standard_deviation),
						sprintf("%3d",$$y->sample_count),
					   ));
	  $self->_print(CRLF);

	  push @outliers, [$j,$i] if $$y->is_outlier;
	  push @modified, [$j,$i] if $$y->is_modified;
	  push @masks,    [$j,$i] if $$y->is_masked;

	}
  }

  $self->_print(CRLF);

  $self->_print("[MASKS]" . CRLF . $array->masks);
  foreach my $mask (@masks){
	next if $mask->[0] == 0 and $mask->[1] == 0;  #WHY ???
	$self->_print($mask->[0], "\t", $mask->[1], CRLF);
  }
  $self->_print(CRLF);

  $self->_print("[OUTLIERS]" . CRLF . $array->outliers);
  foreach my $outlier (@outliers){
	$self->_print($outlier->[0], "\t", $outlier->[1], CRLF);
  }
  $self->_print(CRLF);

  $self->_print("[MODIFIED]" . CRLF . $array->modified);
  foreach my $modified (@modified){
	$self->_print($modified->[0], "\t", $modified->[1], CRLF);
  }
}

1;
