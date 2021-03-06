# $Id$
# BioPerl module for Bio::Expression::Microarray::Feature
#
# Copyright Allen Day <allenday@ucla.edu>, Stanley Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::Affymetrix::Feature - an Affymetrix feature.

=head1 SYNOPSIS

 #get a Feature object somehow...

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
package Bio::Expression::Microarray::Affymetrix::Feature;

use strict;
use base qw(Bio::Expression::Microarray::FeatureI Bio::Root::Root);
use vars qw($DEBUG);

use Class::MakeMethods::Template::Flyweight
  scalar => [qw(
				probe feat expos pos cbase pbase tbase
				atom index codon_index codon regiontype region
				length value standard_deviation sample_count
				display_id
				x y is_match is_masked is_outlier is_modified is_singleton
			   )
			],
  new => 'new',
;

=head2 new

 Title   : new
 Usage   : $ftr = Bio::Expression::Microarray::Affymetrix::Feature->new();
 Function: create a new feature object
 Returns : a Bio::Expression::Microarray::Affymetrix::Feature object
 Args    : none.  all attributes must be set by calling the
           appropriate method
=cut

=head2 get/set methods

The following methods can be used to set or retrieve an attribute
for a feature object.  Call them as in (a) to set an attribute, or
as in (b) to retrieve the value of an attribute:

 (a) $ftr->method('new_value');
 (b) $ftr->method();

Note that no attempt is made to validate the values you store
using an accessor method.

The following methods, along with brief descriptions of their
purpose, are available:

 Method                   Purpose
 ------                   -------
 probe                    ???
 feat                     ???
 expos                    ???
 pos                      ???
 cbase                    ???
 pbase                    ???
 tbase                    ???
 atom                     ???
 index                    ???
 codon_index              ???
 codon                    ???
 regiontype               ???
 region                   ???

=cut


1;
