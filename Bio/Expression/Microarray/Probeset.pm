# $Id$
# BioPerl module for Bio::Expression::Microarray::Probeset
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelsonley <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::Probeset - a set of DNA microarray probes

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

use base qw(Bio::Root::Root);
use vars qw($DEBUG);

=head2 new

 Title   : new
 Usage   : $probeset = Bio::Expression::Microarray::Probeset->new(%args);
 Function: create a new probeset object
 Returns : a Bio::Expression::Microarray::Probeset object
 Args    : an optional hash of parameters to be used in initialization:
           -id    --  the probeset ID
           -type  --  the probeset type

=cut

sub new {
  my($class,@args) = @_;
  my $self = bless {}, $class;
  $self->_initialize(@args);
  return $self;
}

=head2 _initialize

 Title   : _initialize
 Usage   : $probeset->_initialize(@args);
 Function: initialize the probeset object
 Returns : nothing
 Args    : @args

=cut

sub _initialize{
  my ($self,@args) = @_;
  my %param = @args;

  $self->type($param{-type});
  $self->id($param{-id}    );

  $self->SUPER::_initialize(@args);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);
}

=head2 type

 Title   : type
 Usage   : $probeset->type($optional_arg);
 Function: get/set the type of the probeset
 Comments: this is probably going to be a string like
           "quality control", "mismatch blah blah", etc.
 Returns : the probeset type
 Args    : a new value for the probeset type

=cut

sub type {
  my $self = shift;
  $self->{type} = shift if @_;
  return $self->{type};
}

=head2 id

 Title   : id
 Usage   : $probeset->id($optional_arg);
 Function: get/set the id of the probeset
 Returns : the probeset id
 Args    : a new value for the probeset id

=cut

sub id {
  my $self = shift;
  $self->{id} = shift if @_;
  return $self->{id};
}

=head2 add_probe

 Title   : add_probe
 Usage   : $probe_copy = $probeset->add_probe($probe);
 Function: add a probe to the probeset
 Returns : see this_probe()
 Args    : a Bio::Expression::ProbeI compliant object

=cut

sub add_probe {
  my($self,@args) = @_;
  foreach my $probe (@args){
	$self->throw('Probes must be Bio::Expression::ProbeI compliant') unless $probe->isa('Bio::Expression::ProbeI');
    push @{$self->{probes}}, $probe;
  }

  return $self->{probes} ? $self->{probes}->[-1] : undef;
}

=head2 this_probe

 Title   : this_probe
 Usage   : $probe = $probeset->this_probe
 Function: access the last probe added to the probeset
 Returns : the last probe added to the probeset
 Args    : none

=cut

sub this_probe {
  my $self = shift;
  return $self->{probes} ? $self->{probes}->[-1] : undef;
}

=head2 each_probe

 Title   : each_probe
 Usage   : @probes = $probeset->each_probe
 Function: returns a list of Bio::Expression::ProbeI compliant
           objects
 Returns : a list of objects
 Args    : none

=cut

sub each_probe {
  my $self = shift;
  return @{$self->{probes}};
}

=head2 each_probe_value

 Title   : each_probe_value
 Usage   : @probevalues = $probeset->each_probe_value;
 Function: returns an list of values of the probes in the probeset
 Returns : a list of numeric values
 Args    : none

=cut

sub each_probe_value {
  my $self = shift;
  my @values = ();
  push @values, $_->value foreach $self->each_probe;
  return @values;
}

1;
