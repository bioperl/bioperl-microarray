# $Id$
# BioPerl module for Bio::MicroarrayIO::cel
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Microarray::cel - Affymetrix CEL file.

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
package Bio::Expression::Microarray::Affymetrix::Cel;

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use IO::File;
#use PDL;
#use PDL::Matrix;

use base qw(Bio::Root::Root Bio::Root::IO);
use vars qw($DEBUG);

use Class::MethodMaker
#use Class::MakeMethods::Emulator::MethodMaker
  get_set => [qw( cdf )],
  new_with_init => 'new',
;

sub _initialize {
  return shift->init(@_);
}

sub init{
  my ($self,@args) = @_;
  my %args = map {$_=>1} @args;

  $self->SUPER::_initialize(@args);
  $self->_initialize_io(@args);
  my ($usetempfile) = $self->_rearrange([qw(TEMPFILE)],@args);
  defined $usetempfile && $self->use_tempfile($usetempfile);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);

#warn $self->cdf;
  foreach my $arg (keys %args){
    $self->$_($args{$_}) if $self->can($_);
  }
#warn $self->cdf;

#  $self->load_cel;
}

sub matrix {
  my($self,@args) = @_;
  return $self->{matrix} if( defined $self->{matrix} and !@args );

  $self->{matrix} = [] if !defined $self->{matrix};
  $self->{matrix} = mzeroes(@args) if @args and $self->use_pdl and !defined $self->{matrix};

  return $self->{matrix};
}

sub load_cel {
  my($self,@args) = @_;
  my($tfh);

  if( $self->use_tempfile ) {
	$tfh = IO::File->new_tmpfile or $self->throw("Unable to open temp file: $!");
	$tfh->autoflush(1);
  }
  my $mode = undef;
  while( defined( $_ = $self->_readline ) ){
    my($key,$value) = (undef,undef);
    if(my($try) = /^\[(.+)\]/){
      $mode = $try;
      next;
    } else {
      ($key,$value) = $_ =~ /^(.+?)=(.+)$/;
    }

    if($mode eq 'CEL'){
      $self->{$key} = $value if $key;
    }
    elsif($mode eq 'HEADER'){
      $self->{$key} = $value if $key;
    }
    elsif($mode eq 'INTENSITY'){
      if($key){
        $self->{$key} = $value if $key;
        next;
      }

      s/\s*(.+)\s*/$1/; #clean up spaces;
      my @row = split /\s+/;

print STDERR sprintf("%-60s\r",sprintf("%3d %3d is %f",$row[1],$row[0],$row[2]));

	my $probe = $self->cdf->matrix($row[0],$row[1]);
	$$probe->value($row[2]) if $probe;

    }
    elsif($mode eq 'MASKS'){}
    elsif($mode eq 'OUTLIERS'){}
    elsif($mode eq 'MODIFIED'){}
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
