# $Id$
# BioPerl module for Bio::Expression::Microarray::Probeset
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::Probeset - a cluster of DNA microarray probes

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
package Bio::Expression::Microarray::Probeset;

use strict;
#use PDL;
#use PDL::Matrix;

use base qw(Bio::Root::Root);
use vars qw($DEBUG);

sub new {
  my($class,@args) = @_;
  my $self = bless {}, $class;
  $self->_initialize(@args);
  return $self;
}

sub _initialize{
  my ($self,@args) = @_;

  $self->SUPER::_initialize(@args);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);
}

sub type {
  my $self = shift;
  $self->{type} = shift if @_;
  return $self->{type};
}

sub id {
  my $self = shift;
  $self->{id} = shift if @_;
  return $self->{id};
}

sub add_probe {
  my($self,@args) = @_;
  foreach my $probe (@args){
	$self->throw('Probes must be Bio::Expression::ProbeI compliant') unless $probe->isa('Bio::Expression::ProbeI');
    push @{$self->{probes}}, $probe;
  }

  return $self->{probes} ? @{$self->{probes}} : undef;
}

sub this_probe {
  my $self = shift;
  return $self->{probes} ? $self->{probes}->[-1] : undef;
}

sub each_probe {
  my $self = shift;
  return @{$self->{probes}};
}

sub each_probe_value {
  my $self = shift;
  my @values = ();
  push @values, $_->value foreach $self->each_probe;
  return @values;
}

sub each_probe_value_PDL {
  my $self = shift;
  eval{require PDL};
  if($@){ $self->throw("couldn't load PDL: $@") }

  return PDL->new($self->values);
}

1;
