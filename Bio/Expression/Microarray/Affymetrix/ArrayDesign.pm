# $Id$
# BioPerl module for Bio::Expression::Affymetrix::ArrayDesign
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expresssion::Affymetrix::ArrayDesign - Affy Chip Template.

=head1 SYNOPSIS


=head1 DESCRIPTION


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
package Bio::Expression::Microarray::Affymetrix::ArrayDesign;

use strict;
use Bio::Root::Root;
use Bio::Expression::FeatureGroup;
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

=head2 _initialize

 Title   : _initialize
 Function: For compatibility with Bioperl. Defers to init().

=cut

sub _initialize {
  return shift->init(@_);
}

=head2 init

 Title   : init
 Function: For compatibility with Class::MakeMethods. 

=cut

sub init {
  my ($self,@args) = @_;

  $self->SUPER::_initialize(@args);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);
}

=head2 matrix

 Title   : matrix
 Usage   : $self->matrix($x_coord, $y_coord, \$feature)
 Function: get-set method for matrix object/coordinate pair
 Returns : a matrix location
 Args    : A coordinate for a matrix location and optional value

=cut

sub matrix {
  my($self,@args) = @_;
  $self->{matrix} = []   unless defined $self->{matrix};
  return $self->{matrix} unless defined $args[0];

  $self->{matrix}->[ $args[1] ][ $args[0] ] = $args[2] if defined $args[2];
  return $self->{matrix}->[ $args[1] ][ $args[0] ];
}

=head2 featuregroup

 Title   : featuregroup
 Usage   : $self->featuregroup->($featurename);
 Function: get-set method for FeatureGroup object
Returns : a Bio::Expression::FeatureGroup object
 Args    : A key for a FeatureGroup object

=cut

sub featuregroup {
  my($self,$arg) = @_;
  return $self->{featuregroup}->{$arg} if $self->{featuregroup}->{$arg};
  $self->{featuregroup}->{$arg} = Bio::Expression::FeatureGroup->new()
	or $self->throw("Couldn't create a Bio::Expression::FeatureGroup: $!");
  return $self->{featuregroup}->{$arg};
}

=head2 qc_featuregroup

 Title   : qc_featuregroup
 Usage   : $self->qc_featuregroup($mode);
 Function: get-set method for quality control FeatureGroup object
Returns : a Bio::Expression::FeatureGroup object
 Args    : A key for a FeatureGroup object

=cut

sub qc_featuregroup {
  my($self,$arg) = @_;
  return $self->{qcfeaturegroup}->{$arg} if $self->{qcfeaturegroup}->{$arg};
  $self->{qcfeaturegroup}->{$arg} = Bio::Expression::FeatureGroup->new()
	or $self->throw("Couldn't create a Bio::Expression::FeatureGroup: $!");

  #tag it as being a QC featuregroup
  $self->{qcfeaturegroup}->{$arg}->is_qc(1);

  return $self->{qcfeaturegroup}->{$arg};
}

=head2 each_featuregroup

 Title   : each_featuregroup
 Usage   : @featuregroups = $array->each_featuregroup();
 Function: gets a list of FeatureGroup objects
 Returns : returns list of FeatureGroup objects
 Args    : none

=cut

sub each_featuregroup {
  my $self = shift;
  my @return = ();
  foreach my $p (sort keys %{$self->{featuregroup}}){
	push @return, $self->{featuregroup}->{$p};
  }
  return @return;
}

=head2

 Title   : each_qcfeaturegroup
 Usage   : @qcfeaturegroups = $array->each_qcfeaturegroup();
 Function: gets a list of quality control FeatureGroup objects
 Returns : returns list of quality control FeatureGroup objects
 Args    : none

=cut

sub each_qcfeaturegroup {
  my $self = shift;
  my @return = ();
  foreach my $p (sort keys %{$self->{qcfeaturegroup}}){
	push @return, $self->{qcfeaturegroup}->{$p};
  }
  return @return;
}

=head2

 Title   : load_data
 Usage   : $array->load_data($line);
 Function: parses current line of file and loads information
 Returns : nothing
 Args    : The line of text to be parsed

=cut
sub load_data {
  my($self,$line) = @_;

  return unless $line;
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

	my $featuregroup = $self->qc_featuregroup($self->mode);

	my($type) = $_ =~ /Type=(.+)/;

	$featuregroup->type($type) and return if $type;
	$featuregroup->id($self->mode) if $self->mode;

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
	$featuregroup->add_feature($feature);
  }
  elsif($self->mode =~ /^Unit(\d+)_Block/){
	return if /^Block|Num|Start|Stop|CellHeader/;

	my $featuregroup;

	my($name) = $_ =~ /^Name=(.+)/;
	if($name){
	  $featuregroup = $self->featuregroup($name);
	  $featuregroup->id($name);
	  $self->_temp_name($name);
	  return;
	} else {
	  $featuregroup = $self->featuregroup($self->_temp_name);
	}

	my($feature,$attrs) = (undef,undef);
	($feature,$attrs) = $_ =~ /Cell(\d+)=(.+)/;
	my @attrs = split /\t/, $attrs;

	my %featureparams = (
						 x        =>	$attrs[UNIT_X],
						 id       =>	$attrs[UNIT_QUAL],
						 y	    =>	$attrs[UNIT_Y],
						 is_match =>  $attrs[UNIT_CBASE] eq $attrs[UNIT_PBASE] ? 0 : 1,
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
	$featuregroup->add_feature($feature);

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

# class::makemethods::template::flyweight no longer implements destroy_flyweight_info ?
#  foreach my $x (@{$self->matrix}){
#    next unless $x;
#    foreach my $y (@{$x}){
#      next unless $y;
#      $$y->_destroy_flyweight_info;
#    }
#  }
}

1;
