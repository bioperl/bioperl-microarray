# $Id$
# BioPerl module for Bio::Expression::Affymetrix::dChipXLS
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expresssion::Affymetrix::dChipArray - dChip processed Affymetrix data.

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
package Bio::Expression::Microarray::Affymetrix::dChipArray;

use strict;
use Bio::Root::Root;
use Bio::Expression::FeatureGroup;
use Bio::Expression::Microarray::Affymetrix::Feature;

use base qw(Bio::Root::Root);
use vars qw($DEBUG);

use Class::MakeMethods::Emulator::MethodMaker
  get_set       => [ qw( id
					   )
				   ],
  new_with_init => 'new',
;

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: get/set the ID of the feature
 Returns : value of id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

=head2 _initialize

 Title   : _initialize
 Usage   : this is for Bioperl compatibility, use init() instead

=cut

sub _initialize {
  return shift->init(@_);
}

=head2 init

 Title   : init
 Usage   : $obj->init(@args)
 Function:
 Example :
 Returns :
 Args    :


=cut

sub init {
  my ($self,@args) = @_;

  $self->SUPER::_initialize(@args);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);
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

=head2 load_data

 
 Title   : load_data
 Usage   : $dchip->load_data($line);
 Function: parses line from dchip file and loads data
 Returns : nothing
 Args    : The line of the parsed file to be loaded

=cut

sub load_data {
  my($self,$line) = @_;

  next unless $line;
  print STDERR $self->mode . "\r" if $DEBUG;

  my($probe,$value,$call,$sd) = split /\t/, $line;

  if($probe =~ /^probe ?set$/ and !$self->id()){
    $self->id($value);
    return;
  }

  my $featuregroup = $self->featuregroup($probe);
  $featuregroup->id($probe);
  $featuregroup->standard_deviation($sd);
  $featuregroup->quantitation($value);
  $featuregroup->presence($call);
}

1;
