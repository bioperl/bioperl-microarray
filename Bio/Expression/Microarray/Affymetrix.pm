# $Id$
# BioPerl module for Bio::Expression::Microarray::Affymetrix
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::Affymetrix - Work with an Affymetrix array.

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
package Bio::Expression::Microarray::Affymetrix;

use strict;
use Bio::Root::Root;
use Bio::Expression::Microarray::Affymetrix::Cel;
use Bio::Expression::Microarray::Affymetrix::CDF;

use base qw(Bio::Root::Root);
use vars qw($DEBUG);

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
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);

  my %param = @args;
  @param{ map { lc $_ } keys %param } = values %param; # lowercase keys

  $self->celfile($param{-celfile}) || die "no CEL parameter provided to new()";
  $self->cdffile($param{-cdffile}) || die "no CDF parameter provided to new()";

#warn "loading cdf...";
  $self->load_cdf();
#warn "loading cel...";
  $self->load_cel();
}

sub celfile {
  my($self,$val) = @_;
  $self->{celfile} = $val if $val;
  return $self->{celfile};
}

sub cel {
  my($self,$val) = @_;
  $self->{cel} = $val if $val;
  return $self->{cel};
}

sub cdffile {
  my($self,$val) = @_;
  $self->{cdffile} = $val if $val;
  return $self->{cdffile};
}

sub cdf {
  my($self,$val) = @_;
  $self->{cdf} = $val if $val;
  return $self->{cdf};
}

sub load_cel {
  my($self,@args) = @_;
  my($tfh);

  my $cel = Bio::Expression::Microarray::Affymetrix::Cel->new(
						  -file => $self->celfile,
						  );
  $cel->cdf($self->cdf);

  $cel->load_cel;
  $self->cel($cel);
}

sub load_cdf {
  my($self,@args) = @_;
  my($tfh);

  my $cdf = Bio::Expression::Microarray::Affymetrix::CDF->new(
						  -file => $self->cdffile,
#			                           cel  => $self->cel,
						  );
  $cdf->load_cdf;
  $self->cdf($cdf);
}

1;
