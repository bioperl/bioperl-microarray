# $Id$
# BioPerl module for Bio::Expression::Microarray::Probe
#
# Copyright Allen Day <allenday@ucla.edu>, Stanley Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::Affymetrix::Probe - an Affymetrix probe.

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
package Bio::Expression::Microarray::Affymetrix::Probe;

use strict;
use base qw(Bio::Expression::Microarray::ProbeI Bio::Root::Root);
use vars qw($DEBUG);

use Class::MakeMethods::Template::Flyweight
  scalar => [qw(
				probe feat expos pos cbase pbase tbase
				atom index codon_index codon regiontype region
				length value standard_deviation sample_count
			   )
			],
  new => 'new',
;

=head2 new

 Title   : new
 Usage   : $probe = Bio::Expression::Microarray::Affymetrix::Probe->new();
 Function: create a new probe object
 Returns : a Bio::Expression::Microarray::Affymetrix::Probe object
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
