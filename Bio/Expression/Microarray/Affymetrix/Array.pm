# $Id$
# BioPerl module for Bio::Expression::Affymetrix::Array
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expresssion::Affymetrix::Array - Affy Chip Template.

=head1 SYNOPSIS


=head1 DESCRIPTION


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
package Bio::Expression::Microarray::Affymetrix::Array;

use strict;
use Bio::Root::Root;
use Bio::Expression::FeatureSet;
use Bio::Expression::Microarray::Affymetrix::Feature;

use base qw(Bio::Root::Root);
use vars qw($DEBUG);
use enum qw(:QC_ X Y PROBE PLEN ATOM INDEX MATCH BG);
use enum qw(:UNIT_ X Y PROBE FEAT QUAL EXPOS POS CBASE PBASE TBASE ATOM INDEX CODONIND CODON REGIONTYPE REGION);

use Class::MakeMethods::Emulator::MethodMaker
  get_set       => [ qw(
						cel header modified intensity masks outliers heavy
						algorithm algorithm_parameters name date type version dat_header
						mode _temp_name id
					   )
				   ],
  new_with_init => 'new',
;

sub _initialize {
  return shift->init(@_);
}

sub init {
  my ($self,@args) = @_;

  $self->SUPER::_initialize(@args);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);
}

sub matrix {
  my($self,@args) = @_;
  $self->{matrix} = []   unless defined $self->{matrix};
  return $self->{matrix} unless defined $args[0];

  $self->{matrix}->[ $args[1] ][ $args[0] ] = $args[2] if defined $args[2];
  return $self->{matrix}->[ $args[1] ][ $args[0] ];
}

sub featureset {
  my($self,$arg) = @_;
  return $self->{featureset}->{$arg} if $self->{featureset}->{$arg};
  $self->{featureset}->{$arg} = Bio::Expression::FeatureSet->new()
	or $self->throw("Couldn't create a Bio::Expression::FeatureSet: $!");
  return $self->{featureset}->{$arg};
}

sub qcfeatureset {
  my($self,$arg) = @_;
  return $self->{qcfeatureset}->{$arg} if $self->{qcfeatureset}->{$arg};
  $self->{qcfeatureset}->{$arg} = Bio::Expression::FeatureSet->new()
	or $self->throw("Couldn't create a Bio::Expression::FeatureSet: $!");

  #tag it as being a QC featureset
  $self->{qcfeatureset}->{$arg}->is_qc(1);

  return $self->{qcfeatureset}->{$arg};
}

sub each_featureset {
  my $self = shift;
  my @return = ();
  foreach my $p (sort keys %{$self->{featureset}}){
	push @return, $self->{featureset}->{$p};
  }
  return @return;
}

sub each_qcfeatureset {
  my $self = shift;
  my @return = ();
  foreach my $p (sort keys %{$self->{qcfeatureset}}){
	push @return, $self->{qcfeatureset}->{$p};
  }
  return @return;
}

sub load_data {
  my($self,$line) = @_;

  next unless $line;
  print STDERR $self->mode . "\r" if $DEBUG;
  my($key,$value) = (undef,undef);

  if(my($try) = $line =~ /^\[(.+)\]/){
	$self->mode($try);
	return;
  } else {
	($key,$value) = $line =~ /^(.+?)=(.+)$/;
  }

  if($self->mode eq 'CDF'){
	$self->{lc($key)} = $value if $key;
  }
  elsif($self->mode eq 'Chip'){
	$self->{lc($key)} = $value if $key;
  }

  elsif($self->mode =~ /^QC/){
	return if /^CellHeader/;

	my $featureset = $self->qcfeatureset($self->mode);

	my($type) = $_ =~ /Type=(.+)/;

	$featureset->type($type) and return if $type;
	$featureset->id($self->mode) if $self->mode;

	my($feature,$attrs) = $_ =~ /Cell(\d+)=(.+)/;
	return unless $attrs;
	my @attrs = split /\t/, $attrs;

	my %featureparams = (
						 x		=>	$attrs[QC_X],
						 y		=>	$attrs[QC_Y],
						);

	if($self->heavy){
	  $featureparams{probe}  = 	$attrs[QC_PROBE];
	  $featureparams{length} = 	$attrs[QC_PLEN];
	  $featureparams{atom}   = 	$attrs[QC_ATOM];
	  $featureparams{index}  = 	$attrs[QC_INDEX];
	}

	my $feature = Bio::Expression::Microarray::Affymetrix::Feature->new( %featureparams );
	$self->matrix($attrs[UNIT_X],$attrs[UNIT_Y],\$feature);
	$featureset->add_feature($feature);
  }
  elsif($self->mode =~ /^Unit(\d+)_Block/){
	return if /^Block|Num|Start|Stop|CellHeader/;

	my $featureset;

	my($name) = $_ =~ /^Name=(.+)/;
	if($name){
	  $featureset = $self->featureset($name);
	  $featureset->id($name);
	  $self->_temp_name($name);
	  return;
	} else {
	  $featureset = $self->featureset($self->_temp_name);
	}

	my($feature,$attrs) = $_ =~ /Cell(\d+)=(.+)/;
	my @attrs = split /\t/, $attrs;

	my %featureparams = (
						 x        =>	$attrs[UNIT_X],
						 id       =>	$attrs[UNIT_QUAL],
						 y	    =>	$attrs[UNIT_Y],
						 is_match =>  $attrs[UNIT_CBASE] eq $attrs[UNIT_PBASE] ? 1 : 0,
						);

	if($self->heavy){
	  $featureparams{probe}		= 	$attrs[UNIT_PROBE];
	  $featureparams{feat}		= 	$attrs[UNIT_FEAT];
	  $featureparams{expos}		= 	$attrs[UNIT_EXPOS];
	  $featureparams{pos}		= 	$attrs[UNIT_POS];
	  $featureparams{cbase}		= 	$attrs[UNIT_CBASE];
	  $featureparams{pbase}		= 	$attrs[UNIT_PBASE];
	  $featureparams{tbase}		= 	$attrs[UNIT_TBASE];
	  $featureparams{atom}		= 	$attrs[UNIT_ATOM];
	  $featureparams{index}		= 	$attrs[UNIT_INDEX];
	  $featureparams{codon_index}	= 	$attrs[UNIT_CODONIND];
	  $featureparams{codon}		= 	$attrs[UNIT_CODON];
	  $featureparams{regiontype}	= 	$attrs[UNIT_REGIONTYPE];
	  $featureparams{region}		= 	$attrs[UNIT_REGION];
	}

	my $feature = Bio::Expression::Microarray::Affymetrix::Feature->new( %featureparams );
	$featureset->add_feature($feature);

	$self->matrix($attrs[UNIT_X],$attrs[UNIT_Y],\$feature);
    }
  elsif($self->mode =~ /^Unit(\d+)/){
	#not sure what should be done with these... they seem extraneous
  }
}

sub DESTROY {
  my $self = shift;
  $self->destroy_features();
}

sub destroy_features {
  my $self = shift;
  my $matrix = $self->matrix;

  foreach my $x (@{$self->matrix}){
    next unless $x;
    foreach my $y (@{$x}){
      next unless $y;
      $$y->_destroy_flyweight_info;
    }
  }
}

1;
