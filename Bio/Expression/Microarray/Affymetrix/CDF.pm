# $Id$
# BioPerl module for Bio::Expression::Microarray::Affymetrix::CDF
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expresssion::Microarray::Affymetrix::CDF - Affy Chip Definition File.

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
package Bio::Expression::Microarray::Affymetrix::CDF;

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
#use PDL;
#use PDL::Matrix;
use Bio::Expression::Microarray::Probeset;
use Bio::Expression::Microarray::Probe;

use base qw(Bio::Root::Root Bio::Root::IO);
use vars qw($DEBUG);
use enum qw(:QC_ X Y PROBE PLEN ATOM INDEX MATCH BG);
use enum qw(:UNIT_ X Y PROBE FEAT QUAL EXPOS POS CBASE PBASE TBASE ATOM INDEX CODONIND CODON REGIONTYPE REGION);

sub new {
  my($class,@args) = @_;
  my $self = bless {}, $class;
  $self->_initialize(@args);
  return $self;
}

sub _initialize{
  my ($self,@args) = @_;

#  my %param = @args;
#  $self->$_($param{$_}) foreach (keys %param);

  $self->SUPER::_initialize(@args);
  $self->_initialize_io(@args);
  my ($usetempfile) = $self->_rearrange([qw(TEMPFILE)],@args);
  defined $usetempfile && $self->use_tempfile($usetempfile);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);

  $self->load_cdf;
}

sub matrix {
  my($self,@args) = @_;
  $self->{matrix} = [] unless defined $self->{matrix};

  $self->{matrix}->[$args[0]][$args[1]] = $args[2] if defined $args[2];
  return $self->{matrix}->[$args[0]][$args[1]];
}

sub heavy {
  my($self,$arg) = @_;
  return $self->{heavy} if !$arg;
  $self->{heavy} = $arg;
  return $self->{heavy};
}

sub probeset {
  my($self,$arg) = @_;
  return $self->{probeset}->{$arg} if $self->{probeset}->{$arg};
  $self->{probeset}->{$arg} = Bio::Expression::Microarray::Probeset->new()
	or $self->throw("Couldn't create a Bio::Expression::Microarray::Probeset: $!");
  return $self->{probeset}->{$arg};
}

sub qcprobeset {
  my($self,$arg) = @_;
  return $self->{qcprobeset}->{$arg} if $self->{qcprobeset}->{$arg};
  $self->{qcprobeset}->{$arg} = Bio::Expression::Microarray::Probeset->new()
	or $self->throw("Couldn't create a Bio::Expression::Microarray::Probeset: $!");
  return $self->{qcprobeset}->{$arg};
}

sub each_probeset {
  my $self = shift;
  my @return = ();
  foreach my $p (sort keys %{$self->{probeset}}){
	push @return, $self->{probeset}->{$p};
  }
  return @return;
}

sub each_qcprobeset {
  my $self = shift;
  my @return = ();
  foreach my $p (sort keys %{$self->{qcprobeset}}){
	push @return, $self->{qcprobeset}->{$p};
  }
  return @return;
}

sub load_cdf {
  my($self,@args) = @_;
  my($tfh);

  my %unit_name = (); #map unit blocks to probe names

  if( $self->use_tempfile ) {
	$tfh = $self->tempfile() or self->throw("Unable to open tempfile: $!");
	$tfh->autoflush(1);
  }

  my $mode = undef;
  my %current = ();

  while( defined( $_ = $self->_readline ) ){
	chomp;
	next unless $_;

print "$mode\r";

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

      my $probeset = $self->qcprobeset($mode);

      my($type) = $_ =~ /Type=(.+)/;

      $probeset->type($type) and next if $type;
      $probeset->id($mode) if $mode;

      my($probe,$attrs) = $_ =~ /Cell(\d+)=(.+)/;
	  next unless $attrs;
      my @attrs = split /\t/, $attrs;

          my %probeparams = (
		  x		=>	$attrs[QC_X],
		  y		=>	$attrs[QC_Y],
          );

          if($self->heavy){
                  $probeparams{probe} = 	$attrs[QC_PROBE];
                  $probeparams{length} = 	$attrs[QC_PLEN];
                  $probeparams{atom} = 		$attrs[QC_ATOM];
                  $probeparams{index} = 	$attrs[QC_INDEX];
          }

	  my $probe = Bio::Expression::Microarray::Probe->new( %probeparams );
	  $self->matrix($attrs[UNIT_X],$attrs[UNIT_Y],\$probe);
	  $probeset->add_probe($probe);
    }
    elsif($mode =~ /^Unit(\d+)_Block/){
      next if /^Block|Num|Start|Stop|CellHeader/;

	  my $probeset;

	  my($name) = $_ =~ /^Name=(.+)/;
	  if($name){
		$probeset = $self->probeset($name);
		$probeset->id($name);
		$current{name} = $name;
		next;
	  } else {
		$probeset = $self->probeset($current{name});
	  }

      my($probe,$attrs) = $_ =~ /Cell(\d+)=(.+)/;
      my @attrs = split /\t/, $attrs;

      my %probeparams = (
	   x      =>	$attrs[UNIT_X],
	   name   =>    $attrs[UNIT_QUAL],
	   y	  =>	$attrs[UNIT_Y],
      );

      if($self->heavy){
		$probeparams{probe} = 		$attrs[UNIT_PROBE];
		$probeparams{feat} = 		$attrs[UNIT_FEAT];
		$probeparams{expos} = 		$attrs[UNIT_EXPOS];
		$probeparams{pos} = 		$attrs[UNIT_POS];
		$probeparams{cbase} = 		$attrs[UNIT_CBASE];
		$probeparams{pbase} = 		$attrs[UNIT_PBASE];
		$probeparams{tbase} = 		$attrs[UNIT_TBASE];
		$probeparams{atom} = 		$attrs[UNIT_ATOM];
		$probeparams{index} = 		$attrs[UNIT_INDEX];
		$probeparams{codon_index} = 	$attrs[UNIT_CODONIND];
		$probeparams{codon} = 		$attrs[UNIT_CODON];
		$probeparams{regiontype} = 	$attrs[UNIT_REGIONTYPE];
		$probeparams{region} = 		$attrs[UNIT_REGION];
      }

      my $probe = Bio::Expression::Microarray::Probe->new( %probeparams);

	  $self->matrix($attrs[UNIT_X],$attrs[UNIT_Y],\$probe);

	  $probeset->add_probe($probe);
    }
    elsif($mode =~ /^Unit(\d+)/){
      #not sure what should be done with these... seem extraneous
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
