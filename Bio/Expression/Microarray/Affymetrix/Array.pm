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
use Bio::Root::IO;
use Bio::Expression::FeatureSet;
use Bio::Expression::Microarray::Affymetrix::Feature;

use base qw(Bio::Root::Root Bio::Root::IO);
use vars qw($DEBUG);
use enum qw(:QC_ X Y PROBE PLEN ATOM INDEX MATCH BG);
use enum qw(:UNIT_ X Y PROBE FEAT QUAL EXPOS POS CBASE PBASE TBASE ATOM INDEX CODONIND CODON REGIONTYPE REGION);

use Class::MakeMethods::Emulator::MethodMaker
  get_set       => [ qw(
						cel header modified intensity masks outliers modified heavy
						algorithm algorithm_parameters name date type version dat_header
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
  $self->_initialize_io(@args);
  my ($usetempfile) = $self->_rearrange([qw(TEMPFILE)],@args);
  defined $usetempfile && $self->use_tempfile($usetempfile);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);

  $self->load_cdf;
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

sub load_cdf {
  my($self,@args) = @_;
  my($tfh);

  my %unit_name = (); #map unit blocks to feature names

  if( $self->use_tempfile ) {
	$tfh = $self->tempfile() or self->throw("Unable to open tempfile: $!");
	$tfh->autoflush(1);
  }

  my $mode = undef;
  my %current = ();

  while( defined( $_ = $self->_readline ) ){
	chomp;
	next unless $_;

	print STDERR "$mode\r" if $DEBUG;

    my($key,$value) = (undef,undef);
    if(my($try) = /^\[(.+)\]/){
      $mode = $try;
      next;
    } else {
      ($key,$value) = $_ =~ /^(.+?)=(.+)$/;
    }

    if($mode eq 'CDF'){
      $self->{lc($key)} = $value if $key;
    }
    elsif($mode eq 'Chip'){
      $self->{lc($key)} = $value if $key;
    }

    elsif($mode =~ /^QC/){
      next if /^CellHeader/;

      my $featureset = $self->qcfeatureset($mode);

      my($type) = $_ =~ /Type=(.+)/;

      $featureset->type($type) and next if $type;
      $featureset->id($mode) if $mode;

      my($feature,$attrs) = $_ =~ /Cell(\d+)=(.+)/;
	  next unless $attrs;
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
    elsif($mode =~ /^Unit(\d+)_Block/){
      next if /^Block|Num|Start|Stop|CellHeader/;

	  my $featureset;

	  my($name) = $_ =~ /^Name=(.+)/;
	  if($name){
		$featureset = $self->featureset($name);
		$featureset->id($name);
		$current{name} = $name;
		next;
	  } else {
		$featureset = $self->featureset($current{name});
	  }

      my($feature,$attrs) = $_ =~ /Cell(\d+)=(.+)/;
      my @attrs = split /\t/, $attrs;

      my %featureparams = (
	   x        =>	$attrs[UNIT_X],
	   id       =>	$attrs[UNIT_QUAL],
	   y	    =>	$attrs[UNIT_Y],
	   is_match =>  $attrs[UNIT_CBASE] eq $attrs[UNIT_PBASE] ? 1 : 0,
      );

print STDERR $featureparams{x} . "         " . $featureparams{y} . "\n" if $mode eq 'Unit943_Block1';

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
    elsif($mode =~ /^Unit(\d+)/){
      #not sure what should be done with these... they seem extraneous
    }
  }
}

=head2 use_tempfile

 Title   : use_tempfile
 Usage   : $obj->use_tempfile($newval)
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

1;
