# $Id$
# BioPerl module for Bio::MicroarrayIO::cel
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Microarray::Affymetrix::Data - Affymetrix CEL file.

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
package Bio::Expression::Microarray::Affymetrix::Data;

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Expression::Microarray::Affymetrix::Probe;
use IO::File;

use base qw(Bio::Root::Root Bio::Root::IO);
use vars qw($DEBUG);

use Class::MakeMethods::Emulator::MethodMaker
  get_set => [qw( template mode )],
  new     => 'new',
;

sub load_data {
  my($self,$line) = @_;

  my($key,$value) = (undef,undef);
  if(my($try) = $line =~ /^\[(.+)\]/){
	$self->mode($try);
	next;
  }

  ($key,$value) = $line =~ /^(.+?)=(.+)$/;

  if($self->mode eq 'CEL'){
	$self->template->cel($self->template->cel . $line) if $key;
  }
  elsif($self->mode eq 'HEADER'){
	$self->template->header($self->template->header . $line) if $key;
  }
  elsif($self->mode eq 'INTENSITY'){
	if($key){
	  $self->template->intensity($self->template->intensity . $line);
	  next;
	}
	
	$line =~ s/\s*(.+)\s*/$1/; #clean up spaces;
	my @row = split /\s+/, $line;
	
	print STDERR sprintf("%-60s\r",sprintf("%3d %3d is %4.1f",$row[0],$row[1],$row[2])) if $DEBUG;
	
	next unless defined $row[1] and defined $row[2];
	my $probe = $self->template->matrix($row[0],$row[1]);
	
	if(defined $probe){
	  $$probe->value($row[2]);
	  $$probe->standard_deviation($row[3]);
	  $$probe->sample_count($row[4]);
	} else {
	  $probe = Bio::Expression::Microarray::Affymetrix::Probe->new(
													   x =>	$row[0],
													   y =>	$row[1],
													  );
	  $self->template->matrix($row[0],$row[1],\$probe);
	  $probe->value($row[2]);
	  $probe->standard_deviation($row[3]);
	  $probe->sample_count($row[4]);
	}
  }
  elsif($self->mode eq 'MASKS'){
	$self->template->masks($self->template->masks . $line) and next if $key;

	my @row = split /\t/, $line;
	next unless @row;

	my $probe = $self->template->matrix($row[0],$row[1]);

	if(defined $probe){
	  next if $row[0] == 0 and $row[1] == 0; #why do i need to do this?

	  $$probe->is_masked(1);
	} else {
	  $probe = Bio::Expression::Microarray::Affymetrix::Probe->new(
													   x =>	$row[0],
													   y =>	$row[1],
													  );
	  $self->template->matrix($row[0],$row[1],\$probe);
	  $probe->is_masked(1);
	}
  }
  elsif($self->mode eq 'OUTLIERS'){
	$self->template->outliers($self->template->outliers . $line) and next if $key;
	
	my @row = split /\t/, $line;
	
	my $probe = $self->template->matrix($row[0],$row[1]);
	
	if(defined $probe){
	  $$probe->is_outlier(1);
	} else {
	  $probe = Bio::Expression::Microarray::Affymetrix::Probe->new(
													   x =>	$row[0],
													   y =>	$row[1],
													  );
	  $self->template->matrix($row[0],$row[1],\$probe);
	  $probe->is_outlier(1) and next;
	}
  }
  elsif($self->mode eq 'MODIFIED'){
	$self->template->modified($self->template->modified . $line) and next if $key;
	
	my @row = split /\t/, $line;
	my $probe = $self->template->matrix($row[0],$row[1]);
	
	if(defined $probe){
	  $$probe->is_modified(1);
	} else {
	  $probe = Bio::Expression::Microarray::Affymetrix::Probe->new(
													   x =>	$row[0],
													   y =>	$row[1],
													  );
	  $self->template->matrix($row[0],$row[1],\$probe);
	  $probe->is_modified(1) and next;
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
