#Let the code begin...
package Bio::Expression::MicroarrayIO::mas50txt;

use strict;
use Bio::Root::Root;
use Bio::Expression::MicroarrayIO;
use Bio::Expression::Microarray::Affymetrix::Mas50Data;
#use Bio::Expression::Microarray::Affymetrix::Data;
use IO::File;

use base qw(Bio::Root::Root Bio::Expression::MicroarrayIO);
use vars qw($DEBUG $started);

use constant CRLF => "\015\012";

#methods for inserting values
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
#end methods for inserting values 


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
  my ($usetempfile) = $self->_rearrange([qw(TEMPFILE)],@args);
  defined $usetempfile && $self->use_tempfile($usetempfile);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);

  my %param = @args;
  @param{ map { lc $_ } keys %param } = values %param; # lowercase keys

  if($self->mode eq 'r'){
	print STDERR "loading array data...\n" if $self->verbose;

	$self->datafile($param{-file});
  }
}


#this will do one whole array
sub next_array {
  my $self = shift;
  print STDERR "loading data...\n" if $self->verbose;

  #create an object for parsing the array data;
  my $data = new Bio::Expression::Microarray::Affymetrix::Mas50Data;

  while ( defined( $_ = $self->_readline(-raw=>1) ) ) {
	
	if ($_ =~ /Expression/) {
	  $self->mode('HEADERTYPE2');
	  
	  #for sequence of data files
	  if ($started) {
		$self->_pushback($_);
		undef $started;
		return $self->array;
	  }
	  $started = 1;
	  
	  $self->mode('DATATYPE2');
	  while ( defined( $_ = $self->_readline(-raw=>1) ) ) {
		$data->load_data($_);	      
	  }
	} elsif (0) {
	  #add special cases for new TXT file types here
	} else {
	  $self->mode('HEADERTYPE1');

	  #for sequence of data files
	  if ($started) {
		$self->_pushback($_);
		undef $started;
		return $self->array;
	  }
	  $started = 1;

	  $_ = $self->_readline(-raw=>1);
	  chomp;
	  $data->array->featureset($_);
	  $data->array->name($_);

	  while ( defined( $_ = $self->_readline(-raw=>1) ) ) {
		$data->load_data($_);
	  }
	}	
  }

  print STDERR "loaded data!\n" if $self->verbose;

  #and return the array object, complete with loaded data!
  return $data;
}



#verbatim
sub write_array {
  my($self,$array) = @_;
  $self->throw_not_implemented();
}

1;
