# $Id$
# BioPerl module for Bio::Expression::Microarray::Probe
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::Probe - a DNA microarray probe

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
package Bio::Expression::Microarray::Probe;

use strict;
use base qw(Bio::Expression::ProbeI Bio::Root::Root);
use vars qw($DEBUG);

use Class::MakeMethods::Template::Flyweight
  scalar => [qw(x y probe feat name expos pos cbase pbase tbase atom index codon_index codon regiontype region
		length value value_units
				)],
  new => 'new',
;

1;
