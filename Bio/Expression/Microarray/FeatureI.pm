# $Id$
# BioPerl module for Bio::Expression::Microarray::FeatureI
#
# Copyright Allen Day <allenday@ucla.edu>, Stanley Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::FeatureI - an interface class for microarray
features

=head1 SYNOPSIS

...

=head1 DESCRIPTION

Bio::Expression::Microarray::FeatureI is derived from
Bio::Expression::FeatureI, with methods specific to microarray
features.

=head1 FEEDBACK

Direct feedback to E<lt>allenday@ucla.eduE<gt> or to the Bioperl mailing list (see below).

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bioperl.org/wiki/Mailing_lists - About the mailing lists

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
package Bio::Expression::Microarray::FeatureI;

use strict;
use base qw(Bio::Expression::FeatureI Bio::Root::Root);
use vars qw($DEBUG);

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
 x                        the x-axis component of the feature's coordinate
                          position on the array.

 y                        the y-axis component of the feature's coordinate
                          position on the array.

 is_match                 1 if the feature is a reverse-compliment of the
                          sequence it is measuring, 0 or undef if it is not
                          a reverse compliment

 is_masked                1 if the feature has been masked out, i.e. should
                          be excluded from further data analysis, undef or
                          0 if the feature should be used

 is_outlier               1 if the feature's value is dissimilar to the other
                          features in the same of set.  0 or undef if the
                          value is similar to the other features.

 is_modified              1 if the feature's value has been modified after
                          quantitation.  0 or undef if the value has not
                          been modified after quantitation.

=cut

sub x           { shift->throw_not_implemented }
sub y           { shift->throw_not_implemented }
sub is_match    { shift->throw_not_implemented }
sub is_masked   { shift->throw_not_implemented }
sub is_outlier  { shift->throw_not_implemented }
sub is_modified { shift->throw_not_implemented }

1;
