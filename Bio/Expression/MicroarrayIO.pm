# $Id$
#
# BioPerl module for Bio::Expression::MicroarrayIO.pm
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day, Stanley Nelson
# Human Genetics, UCLA Medical School, University of California,
# Los Angeles
#
# You may distribute this module under the same terms as perl
# itself

=head1 NAME

Bio::Expression::MicroarrayIO - Handler for Microarray Formats

=head1 SYNOPSIS

  use Bio::Expression::MicroarrayIO;

  $stream  = Bio::Expression::MicroarrayIO->new(
									  '-file'     => "my.cel",
									  '-template' => "my.cdf",
									  '-format'   => "affymetrix",
											   );

  while ( my $in = $stream->next_array() ) {
	print $in->id() . "\n";
  }

=head1 DESCRIPTION

The Bio::Expression::MicroarrayIO module reads various Microarray
data formats such as Affymetrix CEL and CDF.

=head1 CONSTRUCTORS

=head2 Bio::Expression::MicroarrayIO-E<gt>new()

   $str = Bio::Expression::MicroarrayIO->new(-file     => 'filename',
											 -template => 'template',
											 -format   => $format
											);

The new() class method constructs a new Bio::Expression::MicroarrayIO
object.  The returned object can be used to retrieve or print cluster
objects. new() accepts the following parameters:

=over 4

=item -file

A file path to be opened for reading.

=item -format

Specify the format of the file.  Supported formats include:

   affymetrix		*.cel	Affymetrix CEL files

The format name is case insensitive.  'AFFYMETRIX', 'Affymetrix' and
'affymetrix' are all supported.

=item -template

Affymetrix (and other microarray formats?) files require a template
file that defines the location of probes on an array.  This template
is necessary to match values from a matrix of values from a data file
with sets of probes that are on the array.  Pass the path to the
template file as the -template parameter.

=back

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $array = $str-E<gt>next_array()

Fetch the next array from the stream.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/MailList.shtml      - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Expression::MicroarrayIO;

use strict;

use Bio::Root::Root;
use Bio::Root::IO;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : Bio::Expression::MicroarrayIO->new(-file     => 'path/to/filename',
											  -format   => 'format',
											  -template => 'path/to/template')
 Function: Returns a new microarray stream
 Returns : A Bio::Expression::MicroarrayIO handler.
 Args    : -file     => filename
           -format   => format
           -template => template

=cut

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;

    if( $class =~ /Bio::Expression::MicroarrayIO::(\S+)/ ) {
	  my ($self) = $class->SUPER::new(@args);	
	  $self->_initialize(@args);
	  return $self;
    } else {
	  my %param = @args;
	  @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	  my $format = $param{'-format'} || 
	    $class->_guess_format( $param{-file} || $ARGV[0] );
	  $format = "\L$format";	# normalize capitalization to lower case

	  # normalize capitalization
	  return undef unless( &_load_format_module($format) );
	  return "Bio::Expression::MicroarrayIO::$format"->new(@args);
    }
}


# this is borrowed from SeqIO.
# _initialize is chained for all SeqIO classes

sub _initialize {
    my($self, @args) = @_;
    # initialize the IO part
    $self->_initialize_io(@args);
}

=head2 next_array

 Title   : next_array
 Usage   : $ary = $stream->next_array()
 Function: Reads the next array object from the stream and returns it.
 Returns : a Bio::Expression::MicroarrayI compliant object
 Args    : none


=cut

sub next_array {
   my ($self, $seq) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::Expression::MicroarrayIO object.");
}



# this is borrowed from ClusterIO
=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL MicroarrayIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
  my ($format) = @_;
  my ($module, $load, $m);

  $module = "_<Bio/Expression/MicroarrayIO/$format.pm";
  $load = "Bio/Expression/MicroarrayIO/$format.pm";

  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if ( $@ ) {
    print STDERR "
                  $load: couldn't load $format - for more details on
                  supported formats please see the MicroarrayIO docs
                  Exception $@";
	  return;
  }
  return 1;
}

=head2 _guess_format

 Title   : _guess_format
 Usage   : $stream->_guess_format($filename)
 Function: guess format based on file suffix
 Example :
 Returns : guessed format of filename (lower case)
 Args    : filename

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'affymetrix' if /\.cel$/i;
}

=head2 use_tempfile

 Title   : use_tempfile
 Usage   : $stream->use_tempfile($newval)
 Function: Get/Set boolean flag on whether or not use a tempfile
 Example :
 Returns : value of use_tempfile
 Args    : newvalue (optional)

=cut

sub use_tempfile{
  my ($self,$value) = @_;
  if( defined $value) {
	$self->{'_use_tempfile'} = $value;
  }
  return $self->{'_use_tempfile'};
}

sub DESTROY {
    my $self = shift;
    $self->close();
}



1;
