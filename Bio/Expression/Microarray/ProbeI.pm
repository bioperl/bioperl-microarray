# $Id$
# BioPerl module for Bio::Expression::Microarray::Probe
#
# Copyright Allen Day <allenday@ucla.edu>, Stanley Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::ProbeI - an interface class for microarray
probes, derived from Bio::Expression::ProbeI

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
package Bio::Expression::Microarray::ProbeI;

use strict;
use base qw(Bio::Expression::ProbeI Bio::Root::Root);
use vars qw($DEBUG);

use Class::MakeMethods::Template::Flyweight
  scalar => [qw(x y name is_match is_masked is_outlier is_modified
			   )
			],
  new => 'new',
;

=head2 new

 Title   : new
 Usage   : $probe = Bio::Expression::Microarray::Probe->new();
 Function: create a new probe object
 Returns : a Bio::Expression::Microarray::Probe object
 Args    : none.  all attributes must be set by calling the
           appropriate method
=cut

=head2 get/set methods

The following methods can be used to set or retrieve an attribute
for a probe object.  Call them as in (a) to set an attribute, or
as in (b) to retrieve the value of an attribute:

 (a) $probe->method();
 (b) $probe->method('new_value');

Note that no attempt is made to validate the values you store
using an accessor method.

The following methods, along with brief descriptions of their
purpose, are available:

 Method                   Purpose
 ------                   -------
 x                        the x-axis component of the probe's coordinate
                          position on the array.

 y                        the y-axis component of the probe's coordinate
                          position on the array.

 name                     ???

 is_match                 1 if the probe is a reverse-compliment of the
                          sequence it is measuring, 0 or undef if it is not
                          a reverse compliment

 is_masked                1 if the probe has been masked out, i.e. should
                          be excluded from further data analysis, undef or
                          0 if the probe should be used

 is_outlier               1 if the probe's value is dissimilar to the other
                          probes in the same probeset.  0 or undef if the
                          value is similar to the other probes.

 is_modified              1 if the probe's value has been modified after
                          quantitation.  0 or undef if the value has not
                          been modified after quantitation.


=cut


1;
