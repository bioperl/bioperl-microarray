# $Id$
# BioPerl module for Bio::Expression::MicroarrayIO::affymetrix
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::MicroarrayIO::affymetrix - Work with an Affymetrix array.

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 FEEDBACK

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
use Bio::Expression::Microarray::Affymetrix::Data;
use Bio::Expression::Microarray::Affymetrix::Template;
use IO::File;

use base qw(Bio::Root::Root Bio::Expression::MicroarrayIO);
use vars qw($DEBUG);

use constant CRLF => "\015\012";

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

  $self->templatefile($param{-template}) || die "no template provided to new()";
  $self->datafile($param{-file}) || die "no data file provided to new()";

  $self->load_template();
}

sub templatefile {
  my($self,$val) = @_;
  $self->{templatefile} = $val if $val;
  return $self->{templatefile};
}

sub template {
  my($self,$val) = @_;
  $self->{template} = $val if $val;
  return $self->{template};
}

sub load_template {
  my($self,@args) = @_;

  my $template = Bio::Expression::Microarray::Affymetrix::Template->new(
																		-file => $self->templatefile,
																	   );
  $self->template($template);
}

sub next_array {
  my $self = shift;

  my $array = new Bio::Expression::Microarray::Affymetrix::Data;
  $array->template($self->template);

  my $start = 1;
  while( defined( $_ = $self->_readline ) ){
	$_ =~ s/\n/\r\n/gs; #counteract _readline \r-chopping behavior

	#skip to the beginning of the file
	if($_ =~ m!\[CEL\]!){
	  if($start == 0){
		$array = new Bio::Expression::Microarray::Affymetrix::Data;
		$array->template($self->template);
#		return $self->template;
	  }
	  $start = 0;
	}
	next if $start;

	#slurp up the data
	$array->load_data($_);
  }

  #for the last (or only) file in the input stream
  return $self->template;
}

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

  print "[CEL]"       . CRLF . $array->cel       . CRLF;
  print "[HEADER]"    . CRLF . $array->header    . CRLF;
  print "[INTENSITY]" . CRLF . $array->intensity;

  my $i = 0;
  foreach my $x ( @{ $array->matrix } ){
	next unless $x;
	my $j = 0;
	foreach my $y ( @$x ){
	  next unless $y;
	  #$y is a probe object
	  print join "\t", (sprintf("%3d",$j),
						sprintf("%3d",$i),
						$$y->value,
						$$y->standard_deviation,
						sprintf("%3d",$$y->samples),
					   );
	  print CRLF;

	  push @outliers, [$j,$i] if $$y->is_outlier;
	  push @modified, [$j,$i] if $$y->is_modified;
	  push @masks,    [$j,$i] if $$y->is_masked;

	  $j++;
	}
	$i++;
  }

  print CRLF;

  print "[MASKS]" . CRLF . $array->masks;
  foreach my $mask (@masks){
	print $mask->[0], "\t", $mask->[1], CRLF;
  }
  print CRLF;

  print "[OUTLIERS]" . CRLF . $array->outliers;
  foreach my $outlier (@outliers){
	print $outlier->[0], "\t", $outlier->[1], CRLF;
  }
  print CRLF;

  print "[MODIFIED]" . CRLF . $array->modified;
  foreach my $modified (@modified){
	print $modified->[0], "\t", $modified->[1], CRLF;
  }
  #print CRLF;
}

sub datafile {
  my($self,$val) = @_;
  $self->{datafile} = $val if $val;
  return $self->{datafile};
}

sub array {
  my($self,$val) = @_;
  $self->{array} = $val if $val;
  return $self->{array};
}

1;
