# BioPerl module for Bio::Expression::Affymetrix::mas50txt
#
# Copyright Allen Day <allenday@ucla.edu>, Chris To <crsto@ucla.edu> Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

# Let the code begin...

=head1 NAME

Bio::Expresssion::Affymetrix::mas50txt - Work with MAS 5.0 TXT file.

=head1 SYNOPSIS

You should not be using this module directly.  Try
Bio::Expression::MicroarrayIO instead.

=head1 DESCRIPTION

Bio::Expresssion::Affymetrix::mas50txt parses MAS 5.0 package TXT files.

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
Chris To E<lt>crsto@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#Let the code begin...
package Bio::Expression::MicroarrayIO::mas50txt;

use strict;
use Bio::Root::Root;
use Bio::Expression::MicroarrayIO;
use Bio::Expression::Microarray::Affymetrix::Mas50TxtArray;
use IO::File;

use base qw(Bio::Root::Root Bio::Expression::MicroarrayIO);
use vars qw($started);

use constant CRLF => "\015\012";


sub array {
  my($self,$val) = @_;
  $self->{array} = $val if $val;
  return $self->{array};
}


sub datafile {
  my($self,$val) = @_;
  $self->{datafile} = $val if $val;
  return $self->{datafile};
}
 

=head2 new

 Title   : new
 Usage   : Bio::Expression::MicroarrayIO::mas50txt->new(
											  -file     => 'path/to/filename')
 Comments: You should probably not be instantiating this module directly.
           Use MicroarrayIO instead.
 Args    : -file     => filename


=cut
sub new {
  my ($caller,@args) = @_;
  my $class = ref($caller) || $caller;

  my $self = bless {}, $class;

  $self->_initialize(@args);
  return $self;
}

sub _initialize{
  my ($self,@args) = @_;
  $self->SUPER::_initialize(@args);
  $self->_initialize_io(@args);

  my %param = @args;
  @param{ map { lc $_ } keys %param } = values %param; # lowercase keys

}


=head2 next_array

 Title   : next_array
  Usage   : Bio::Expression::Microarray::Affymetrix::Mas50TxtArray->next_array()
 Example : $affymas5->next_array();
 Function: reads an affy Mas 5.0  data record
 Returns : a  Bio::Expression::Microarray::Affymetrix::Mas50TxtArray object
 Args    : none

=cut

sub next_array {
  my $self = shift;
  print STDERR "loading data...\n" if $self->verbose;

  #create an object for parsing the array data;
  my $data = Bio::Expression::Microarray::Affymetrix::Mas50TxtArray->new(); 
  $self->array($data);
  
  while ( defined( $_ = $self->_readline(-raw=>1) ) ) {
	if($_ =~ '==>') {
	  print STDERR "loaded data!\n" if $self->verbose;
	  next;
	}

	if ($_ =~ /Expression Analysis/) {
	  #$data->mode('HEADERTYPE2');
	  
	  #for sequence of data files
	  if ($started) {
		$self->_pushback($_);
		undef $started;
		#return $self->array;
	  }
	  $started = 1;
	  
	  $data->mode('DATATYPE2');
	  while ( defined( $_ = $self->_readline(-raw=>1) ) ) {
		if($_ =~ '==>') {
		  print STDERR "loaded data!\n" if $self->verbose;
		  last;
		}
		$data->load_data($_);	      
	  }

	  return $self->array;

	} elsif (0) {
	  #add special cases for new TXT file types here
	} else {
	  $data->mode('HEADERTYPE1');

	  chomp;
	  $data->id($_);

	  #for sequence of data files
	  if ($started) {
		$self->_pushback($_);
		undef $started;
		#return $self->array;
	  }
	  $started = 1;

	  while ( defined( $_ = $self->_readline(-raw=>1) ) ) {
		if($_ =~ '==>') {
		  print STDERR "new MAS5\n" if $self->verbose;
		  last;
		}
		$data->load_data($_);
	  }

	  return $self->array;
	}	
  }

  print STDERR "loaded data!\n" if $self->verbose;

  #and return the data object, complete with loaded data!
  #return $data;
}


=head2 write_array

 Title   : write_array
 Usage   : not implemented
 Function: 
 Returns : 
 Args    : 

=cut

sub write_array {
  my($self,$array) = @_;
  $self->throw_not_implemented();
}

1;
