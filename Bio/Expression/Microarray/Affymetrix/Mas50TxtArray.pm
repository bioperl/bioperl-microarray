# BioPerl module for Bio::Expression::Affymetrix::Mas50TxtArray
#
# Copyright Allen Day <allenday@ucla.edu>, Chris To <crsto@ucla.edu> Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

# Let the code begin...

=head1 NAME

Bio::Expresssion::Affymetrix::Mas50TxtArray - Affymetrix MAS 5.0 TXT file.

=head1 SYNOPSIS

You should not be using this module directly.  Try
Bio::Expression::MicroarrayIO instead.

=head1 DESCRIPTION


=head1 FEEDBACK

Direct feedback to E<lt>allenday@ucla.eduE<gt> or to the Bioperl mailing list (see below).

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
Chris To E<lt>crsto@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#Let the code begin...
package Bio::Expression::Microarray::Affymetrix::Mas50TxtArray;

use lib '/home2/crsto/cvsroot/bioperl-live/';

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Expression::FeatureGroup::FeatureGroupMas50;
use IO::File;

use base qw(Bio::Root::Root Bio::Root::IO);
use vars qw($DEBUG);

use Class::MakeMethods::Emulator::MethodMaker
  get_set => [ qw( array mode name id

				   #header info
				   arraylot operatorname sampletype
				   sampledescription project comments
				   reagents reagentlot algorithm 
				   cornerplusavg cornerpluscount
				   cornerminusavg cornerminuscount
				   statisticsone statisticstwo

				   bf alpha1 alpha2 tau gamma1h gamma1l
				   gamma2h gamma2l perturbation tgt nf sf
				   sfgene background noise rawq
				 )],
  new     => 'new',
  ;

=head2 featuregroup

 Title   : featuregroup
 Usage   : $array->featuregroup($name);
 Function: get-set method for FeatureGroup object
Returns : a Bio::Expression::FeatureGroup object
 Args    : A key for a FeatureGroup object

=cut

sub featuregroup {
  my($self,$arg) = @_;
  return $self->{featuregroup}->{$arg} if $self->{featuregroup}->{$arg};
  $self->{featuregroup}->{$arg} = Bio::Expression::FeatureGroup::FeatureGroupMas50->new()
	or $self->throw("Couldn't create a Bio::Expression::FeatureGroup::FeatureGroupMas50 $!");
  return $self->{featuregroup}->{$arg};
}

=head2 each_featuregroup

 Title   : each_featuregroup
 Usage   : @featuregroups = $array->each_featuregroup();
 Function: iterated through FeatureGroup objects
 Returns : a list of Bio::Expression::FeatureGroup::FeatureGroupMas50 object
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


=head2 each_qcfeaturegroup

 Title   : each_qcfeaturegroup
 Caveat  : this method cannot be called for a MAS 5.0 TXT file, as the information
           about which probesets are QC and which are not is not present.

=cut

sub each_qcfeaturegroup {
  shift->throw("each_qcfeaturegroup cannot be called for a MAS 5.0 TXT file");  
}

=head2 load_data

 Title   : load_data
 Usage   : Bio::Expression::Microarray::Affymetrix::Mas50TxtArray->load_data(
																		-line => $line);
 Function: parses a line of text from .txt file and loads the data 
 Example : $data->load_data($line);
 Returns : nothing
 Args    : -line => line of data that will be parsed


=cut

sub load_data {
  my($self,$line) = @_;
  chomp($line);
  
  #for the headers
  if ($self->mode eq 'HEADERTYPE2') {
	my(@value) = split /\t/, $line;
	if (@value[0] =~ 'Analysis') {
	  $self->mode('DATATYPE2');
	}
	return;
      
  } elsif ($self->mode eq 'HEADERTYPE1') {
      
	#need to handle the actual array name
	#along with all other stupdi special cases
	chomp($line);
	$line =~ s/\x0d//g;
	my($key, @value) = split /:/, $line;

	substr(@value[0], 0, 2) =~ s/ //;
	if (0 < length(@value[0])) {

	  if ($key eq 'Probe Array Lot') {
		$self->arraylot(@value[0]);
		return;
	  } elsif ($key eq 'Operator Name') {
		$self->operatorname(@value[0]);
		return;
	  } elsif ($key eq 'Sample Type') {
		$self->sampletype(@value[0]);
		return;
	  } elsif ($key eq 'Sample Descripton') {
		$self->sampledescription(@value[0]);
		return;
	  } elsif ($key eq 'Project') {
		$self->project(@value[0]);
		return;
	  } elsif ($key eq 'Comments') {
		$self->comments(@value[0]);
		return;
	  } elsif ($key eq 'Reagents') {
		$self->reagents(@value[0]);
		return;
	  } elsif ($key eq 'Reagent Lot') {
		$self->reagentlot(@value[0]);
		return;
	  } elsif ($key eq 'Algorithm') {
		$self->algorithm(@value[0]);
		return;
	  } elsif ($key eq 'Corner+ Avg') {
		my @temp = split /, /, @value[0];
		$self->cornerplusavg(@temp[0]);
		$self->cornerpluscount(@value[1]);
		return;
	  } elsif ($key eq 'Corner- Avg') { 
		my @temp = split /, /, @value[0];
		$self->cornerminusavg(@temp[0]);
		$self->cornerminuscount(@value[1]);
		return;
	  }
	}

	@value = split/ /, $line;
      
	if (@value) {
	  if (@value[0] =~ 'Probe') {
		$self->mode('DATATYPE1');
		return;
	  }

	  foreach my $x (@value) {
		my($left, $right) = split /=/, $x;
		if (0 < length($right)) {
		  $left = lc($left);
		  $self->$left($right);
		}
	  }
	  return;
	}


	@value = split/\n/;
	unless (defined(@value)){
	  $self->mode('DATATYPE1');
	}

	return;
      
  } elsif ($self->mode eq 'DATATYPE2') {
	chomp($line);
	my @entries = split /\t/, $line;
	
	if (0 == length(@entries[0])) {
	  return;
	}
	elsif (@entries[0] =~ 'Analysis') {
	  return;
	}
	else {
	  if (@entries[0] == '1') {
		$self->id(@entries[1]);
	  }
	  
	  #max column number is 38 (so 0-37)
	  if (0 < length(@entries[1])) {
		my $featuregroup = $self->featuregroup(@entries[2]);
	     
		if (0 < length(@entries[2])) {
		  $featuregroup->id(@entries[2]);
		}
		if (0 < length(@entries[3])) {
		  $featuregroup->stat_pairs(@entries[3]);
		}
		if (0 < length(@entries[4])) {
		  $featuregroup->stat_pairs_used(@entries[4]);
		}
		if (0 < length(@entries[5])) {
		  $featuregroup->signal(@entries[5]);
		}
		if (0 < length(@entries[6])) {
		  $featuregroup->detection(@entries[6]);
		}
		if (0 < length(@entries[7])) {
		  $featuregroup->detection_p_value(@entries[7]);
		}
		if (0 < length(@entries[8])) {
		  $featuregroup->stat_common_pairs(@entries[8]);
		}
		if (0 < length(@entries[9])) {
		  $featuregroup->signal_log_ratio(@entries[9]);
		}
		if (0 < length(@entries[10])) {
		  $featuregroup->signal_log_ratio_low(@entries[10]);
		}
		if (0 < length(@entries[11])) {
		  $featuregroup->signal_log_ratio_high(@entries[11]);
		}
		if (0 < length(@entries[12])) {
		  $featuregroup->change(@entries[12]);
		}
		if (0 < length(@entries[13])) {
		  $featuregroup->change_p_value(@entries[13]);
		}
		if (0 < length(@entries[14])) {
		  $featuregroup->positive(@entries[14]);
		}
		if (0 < length(@entries[15])) {
		  $featuregroup->negative(@entries[15]);
		}
		if (0 < length(@entries[16])) {
		  $featuregroup->pairs(@entries[16]);
		}
		if (0 < length(@entries[17])) {
		  $featuregroup->pairs_used(@entries[17]);
		}
		if (0 < length(@entries[18])) {
		  $featuregroup->pairs_inavg(@entries[18]);
		}
		if (0 < length(@entries[19])) {
		  $featuregroup->pos_fraction(@entries[19]);
		}
		if (0 < length(@entries[20])) {
		  $featuregroup->log_avg(@entries[20]);
		}
		if (0 < length(@entries[21])) {
		  $featuregroup->pos_neg(@entries[21]);
		}
		if (0 < length(@entries[22])) {
		  $featuregroup->avg_diff(@entries[22]);
		}
		if (0 < length(@entries[23])) {
		  $featuregroup->abs_call(@entries[23]);
		}
		if (0 < length(@entries[24])) {
		  $featuregroup->inc(@entries[24]);
		}
		if (0 < length(@entries[25])) {
		  $featuregroup->dec(@entries[25]);
		}
		if (0 < length(@entries[26])) {
		  $featuregroup->inc_ratio(@entries[26]);
		}
		if (0 < length(@entries[27])) {
		  $featuregroup->dec_ratio(@entries[27]);
		}
		if (0 < length(@entries[28])) {
		  $featuregroup->pos_change(@entries[28]);
		}
		if (0 < length(@entries[29])) {
		  $featuregroup->neg_change(@entries[29]);
		}
		if (0 < length(@entries[30])) {
		  $featuregroup->inc_dec(@entries[30]);
		}
		if (0 < length(@entries[31])) {
		  $featuregroup->dposdneg_ratio(@entries[31]);
		}
		if (0 < length(@entries[32])) {
		  $featuregroup->log_avg_ratio_change(@entries[32]);
		}
		if (0 < length(@entries[33])) {
		  $featuregroup->diff_call(@entries[33]);
		}
		if (0 < length(@entries[34])) {
		  $featuregroup->avg_diff_change(@entries[34]);
		}
		if (0 < length(@entries[35])) {
		  $featuregroup->b_a(@entries[35]);
		}
		if (0 < length(@entries[36])) {
		  $featuregroup->fold_change(@entries[36]);
		}
		if (0 < length(@entries[37])) {
		  $featuregroup->sort_score(@entries[37]);
		}
	  }
	}
  } elsif ($self->mode eq 'DATATYPE1') {
	chomp($line);
	my @entries = split /\t/, $line;
     
	for (my $i = 0;  $i < scalar(@entries); $i++) {
	  chomp(@entries[$i]);
	}

	if (!defined(@entries[0])) {
	  return;
	} elsif (@entries[0] eq 'Probe Set Name') {
	  return;
	} else {
	  if (0 < length(@entries[0])) {
		my $featuregroup = $self->featuregroup(@entries[0]);
 
		if (0 < length(@entries[0])) {
		  $featuregroup->id(@entries[0]);
		}
		if (0 < length(@entries[1])) {
		  $featuregroup->stat_pairs(@entries[1]);
		}
		if (0 < length(@entries[2])) {
		  $featuregroup->stat_pairs_used(@entries[2]);
		}
		if (0 < length(@entries[3])) {
		  $featuregroup->signal(@entries[3]);
		}
		if (0 < length(@entries[4])) {
		  $featuregroup->detection(@entries[4]);
		}
		if (0 < length(@entries[5])) {
		  $featuregroup->detection_p_value(@entries[5]);
		}
		if (0 < length(@entries[6])) {
		  $featuregroup->stat_common_pairs(@entries[6]);
		}
		if (0 < length(@entries[7])) {
		  $featuregroup->signal_log_ratio(@entries[7]);
		}
		if (0 < length(@entries[8])) {
		  $featuregroup->signal_log_ratio_low(@entries[8]);
		}
		if (0 < length(@entries[9])) {
		  $featuregroup->signal_log_ratio_high(@entries[9]);
		}
		if (0 < length(@entries[10])) {
		  $featuregroup->change(@entries[10]);
		}
		if (0 < length(@entries[11])) {
		  $featuregroup->change_p_value(@entries[11]);
		}
		if (0 < length(@entries[12])) {
		  $featuregroup->positive(@entries[12]);
		}
		if (0 < length(@entries[13])) {
		  $featuregroup->negative(@entries[13]);
		}
		if (0 < length(@entries[14])) {
		  $featuregroup->pairs(@entries[14]);
		}
		if (0 < length(@entries[15])) {
		  $featuregroup->pairs_used(@entries[15]);
		}
		if (0 < length(@entries[16])) {
		  $featuregroup->pairs_inavg(@entries[16]);
		}
		if (0 < length(@entries[17])) {
		  $featuregroup->pos_fraction(@entries[17]);
		}
		if (0 < length(@entries[18])) {
		  $featuregroup->log_avg(@entries[18]);
		}
		if (0 < length(@entries[19])) {
		  $featuregroup->pos_neg(@entries[19]);
		}
		if (0 < length(@entries[20])) {
		  $featuregroup->avg_diff(@entries[20]);
		}
		if (0 < length(@entries[21])) {
		  $featuregroup->abs_call(@entries[21]);
		}
		if (0 < length(@entries[22])) {
		  $featuregroup->inc(@entries[22]);
		}
		if (0 < length(@entries[23])) {
		  $featuregroup->dec(@entries[23]);
		}
		if (0 < length(@entries[24])) {
		  $featuregroup->inc_ratio(@entries[24]);
		}
		if (0 < length(@entries[25])) {
		  $featuregroup->dec_ratio(@entries[25]);
		}
		if (0 < length(@entries[26])) {
		  $featuregroup->pos_change(@entries[26]);
		}
		if (0 < length(@entries[27])) {
		  $featuregroup->neg_change(@entries[27]);
		}
		if (0 < length(@entries[28])) {
		  $featuregroup->inc_dec(@entries[28]);
		}
		if (0 < length(@entries[29])) {
		  $featuregroup->dposdneg_ratio(@entries[29]);
		}
		if (0 < length(@entries[30])) {
		  $featuregroup->log_avg_ratio_change(@entries[30]);
		}
		if (0 < length(@entries[31])) {
		  $featuregroup->diff_call(@entries[31]);
		}
		if (0 < length(@entries[32])) {
		  $featuregroup->avg_diff_change(@entries[32]);
		}
		if (0 < length(@entries[33])) {
		  $featuregroup->b_a(@entries[33]);
		}
		if (0 < length(@entries[34])) {
		  $featuregroup->fold_change(@entries[34]);
		}
		if (0 < length(@entries[35])) {
		  $featuregroup->sort_score(@entries[35]);
		}
	  }
	}
     
  }
}
    
1;
