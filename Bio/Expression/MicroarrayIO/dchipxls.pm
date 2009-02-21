# $Id$
# BioPerl module for Bio::Expression::MicroarrayIO::dchipxls
#
# Copyright Allen Day <allenday@ucla.edu>, Stanley Nelson
# <snelson@ucla.edu>.  Human Genetics, UCLA Medical School,
# University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::MicroarrayIO::dchipxls - Work with dChip analysis
of an Affymetrix array.

=head1 SYNOPSIS

  $stream  = Bio::Expression::MicroarrayIO->new(
		'-file'     => "my.xls",
		'-format'   => "dchipxls",
											   );

=head1 DESCRIPTION

Bio::Expression::MicroarrayIO::dchipxls parses dChip XLS files.  These
are a processed (normalized) form of Affymetrix CEL files using the
algorithm described in: I<Li C, Hung Wong W. Model-based analysis of
oligonucleotide arrays: model validation, design issues and standard
error application. Genome Biol. 2001;2(8)>

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
package Bio::Expression::MicroarrayIO::dchipxls;

use strict;
use Bio::Root::Root;
use Bio::Expression::MicroarrayIO;
use Bio::Expression::Microarray::Affymetrix::dChipArray;
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

 Title   : new
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

  if($param{-file} or $param{-fh}){
#  if($self->mode eq 'r'){
	print STDERR "loading array template and data...\n" if $self->verbose;

	$self->datafile($param{-file});
  }
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

=head2 next_array

 Title   : next_array
 Usage   : $affy->next_array();
 Function: reads an Affymetrix data record from $affy->datafile
 Returns : a  Bio::Expression::Microarray::Affymetrix::Template object
           that has been filled with probe values.
 Args    :

=cut

sub next_array {
  my $self = shift;
  print STDERR "loading data...\n" if $self->verbose;

  #create an object for parsing the array data;
  my $array = Bio::Expression::Microarray::Affymetrix::dChipArray->new;

  while ( defined( $_ = $self->_readline(-raw=>1) ) ) {
	#if we see the beginning of a new file
	if ($_ =~ m!^probe ?set!) {
	  print STDERR "new CEL\n" if $self->verbose;
	  my($junk,$name) = split /\t/;

	  #if we have seen another CEL before, _pushback() the beginning
	  #of this new CEL, for the next call to next_array() to handle,
	  #and return the preceding CEL.
	  if($started){
		$self->_pushback($_);
		undef $started;
		return $array;
	  }
	  $array->id($name);
	  $started = 1;
	  next;
	}
	#slurp up the data, line by line
	$array->load_data($_);
  }

  print STDERR "loaded data!\n" if $self->verbose;

  #and return the array object, complete with loaded data!
  return $array;
}

=head2 write_array

 Title   : write_array
 Usage   : $affy->write_array($array);
 Function: write an Affymetrix data record using $array
 Returns : nothing.  prints a lot of text.
 Args    : A Bio::Expression::MicroarrayI compliant object.

=cut

sub write_array {
  my($self,$array) = @_;
  $self->throw_not_implemented();
}

1;
