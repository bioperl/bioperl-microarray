# Let the code begin...
package Bio::Expression::Microarray::Affymetrix::Mas50Data;

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Expression::FeatureSet;
use IO::File;

use base qw(Bio::Root::Root Bio::Root::IO);
use vars qw($DEBUG);

use Class::MakeMethods::Emulator::MethodMaker
  get_set => [ qw( array mode name

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
				   #header info
#  new     => 'new',
  ;


sub new {
  my($class, %arg) = @_;
  print join ' ', keys %arg;
  return bless {}, $class;
}

sub featureset {
  my($self,$arg) = @_;
  return $self->{featureset}->{$arg} if $self->{featureset}->{$arg};
  $self->{featureset}->{$arg} = Bio::Expression::FeatureSet::FeatureSetMas50->new()
	or $self->throw("Couldn't create a Bio::Expression::FeatureSet::FeatureSetMas50 $!");
  return $self->{featureset}->{$arg};
}

sub each_featureset {
  my $self = shift;
  my @return = ();
  foreach my $p (sort keys %{$self->{featureset}}) {
	push @return, $self->{featureset}->{$p};
  }
  return @return;
}

sub load_data {
  my($self,$line) = @_;
  chomp($line);

  #for the headers
  if ($self->mode eq 'HEADERTYPE2') {
  } elsif ($self->mode eq 'HEADERTYPE1') {
	#need to handle the actual array name
	#along with all other stupdi special cases
	chomp($line);
	my($key, @value) = split /:/, $line;
	
	if (@value[0]) {
	  substr(@value[0], 0, 2) =~ s/ //;
	  
	  if ($key eq 'Probe Array Lot') {
		$self->array->arraylot(@value[0]);
		return;
	  } elsif ($key eq 'Operator Name') {
		$self->array->operatorname(@value[0]);
		return;
	  } elsif ($key eq 'Sample Type') {
		$self->array->sampletype(@value[0]);
		return;
	  } elsif ($key eq 'Sample Descripton') {
		$self->array->sampledescription(@value[0]);
		return;
	  } elsif ($key eq 'Project') {
		$self->array->project(@value[0]);
		return;
	  } elsif ($key eq 'Comments') {
		$self->array->comments(@value[0]);
		return;
	  } elsif ($key eq 'Reagents') {
		$self->array->reagents(@value[0]);
		return;
	  } elsif ($key eq 'Reagent Lot') {
		$self->array->reagentlot(@value[0]);
		return;
	  } elsif ($key eq 'Algorithm') {
		$self->array->algorithm(@value[0]);
		return;
	  } elsif ($key eq 'Corner+') {
		my @temp = split /, /, @value[0];
		$self->array->cornerplusavg(@temp[0]);
		$self->array->cornerpluscount(@value[1]);
		return;
	  } elsif ($key eq 'Corner-') {
		my @temp = split /, /, @value[0];
		$self->array->cornerplusavg(@temp[0]);
		$self->array->cornerpluscount(@value[1]);
		return;
	  }
	}
	
	@value = split/ /, $line;
	
	if (@value) {
	  foreach my $x (@value) {
		my($left, $right) = split /=/, $line;
		$self->array->{lc($left)} = $right if(defined($right))
	  }
	  return;
	}

	@value = split/\n/;
	$self->mode('DATATYPE1') if(!defined(@value));

	return;
	
  } elsif ($self->mode eq 'DATATYPE2') {
	my @entries = split /\t/, $line;
     
	if (!defined(@entries[0])) {
	  return;
	} elsif (@entries[0] eq 'Analysis Name') {
	  return;
	} else {
	  #max column number is 38 (so 0-37)
	  if (defined(@entries[1])) {
		my $featureset = $self->array->featureset(@entries[1]);

		if (defined(@entries[2])) {
		  $featureset->probe_set_name(@entries[2]);
		}
		if (defined(@entries[3])) {
		  $featureset->stat_pairs(@entries[3]);
		}
		if (defined(@entries[4])) {
		  $featureset->stat_pairs_used(@entries[4]);
		}
		if (defined(@entries[5])) {
		  $featureset->signal(@entries[5]);
		}
		if (defined(@entries[6])) {
		  $featureset->detection(@entries[6]);
		}
		if (defined(@entries[7])) {
		  $featureset->detection_p-value(@entries[7]);
		}
		if (defined(@entries[8])) {
		  $featureset->stat_common_pairs(@entries[8]);
		}
		if (defined(@entries[9])) {
		  $featureset->signal_log_ratio(@entries[9]);
		}
		if (defined(@entries[10])) {
		  $featureset->signal_log_ratio_low(@entries[10]);
		}
		if (defined(@entries[11])) {
		  $featureset->signal_log_ratio_high(@entries[11]);
		}
		if (defined(@entries[12])) {
		  $featureset->change(@entries[12]);
		}
		if (defined(@entries[13])) {
		  $featureset->change_p-value(@entries[13]);
		}
		if (defined(@entries[14])) {
		  $featureset->positive(@entries[14]);
		}
		if (defined(@entries[15])) {
		  $featureset->negative(@entries[15]);
		}
		if (defined(@entries[16])) {
		  $featureset->pairs(@entries[16]);
		}
		if (defined(@entries[17])) {
		  $featureset->pairs_used(@entries[17]);
		}
		if (defined(@entries[18])) {
		  $featureset->pairs_inavg(@entries[18]);
		}
		if (defined(@entries[19])) {
		  $featureset->pos_fraction(@entries[19]);
		}
		if (defined(@entries[20])) {
		  $featureset->log_avg(@entries[20]);
		}
		if (defined(@entries[21])) {
		  $featureset->pos_neg(@entries[21]);
		}
		if (defined(@entries[22])) {
		  $featureset->avg_diff(@entries[22]);
		}
		if (defined(@entries[23])) {
		  $featureset->abs_call(@entries[23]);
		}
		if (defined(@entries[24])) {
		  $featureset->inc(@entries[24]);
		}
		if (defined(@entries[25])) {
		  $featureset->dec(@entries[25]);
		}
		if (defined(@entries[26])) {
		  $featureset->inc_ratio(@entries[26]);
		}
		if (defined(@entries[27])) {
		  $featureset->dec_ratio(@entries[27]);
		}
		if (defined(@entries[28])) {
		  $featureset->pos_change(@entries[28]);
		}
		if (defined(@entries[29])) {
		  $featureset->neg_change(@entries[29]);
		}
		if (defined(@entries[30])) {
		  $featureset->inc_dec(@entries[30]);
		}
		if (defined(@entries[31])) {
		  $featureset->dpos-dneg_ratio(@entries[31]);
		}
		if (defined(@entries[32])) {
		  $featureset->log_avg_ratio_change(@entries[32]);
		}
		if (defined(@entries[33])) {
		  $featureset->diff_call(@entries[33]);
		}
		if (defined(@entries[34])) {
		  $featureset->avg_diff_change(@entries[34]);
		}
		if (defined(@entries[35])) {
		  $featureset->b_a(@entries[35]);
		}
		if (defined(@entries[36])) {
		  $featureset->fold_change(@entries[36]);
		}
		if (defined(@entries[37])) {
		  $featureset->sort_score(@entries[37]);
		}
	  } elsif ($self->mode eq 'DATATYPE1') {
		chomp($line);
		my @entries = split /\t/, $line;

		if (!defined(@entries[0])) {
		  return;
		} elsif (@entries[0] eq 'Probe Set Name') {
		  return;
		} else {
		  if (defined(@entries[0])) {
			my $featureset = $self->array->featureset($self->array->name);

			if (defined(@entries[0])) { 
			  $featureset->probe_set_name(@entries[0]);
			}
			if (defined(@entries[1])) { 
			  $featureset->stat_pairs(@entries[1]);
			}
			if (defined(@entries[2])) { 
			  $featureset->stat_pairs_used(@entries[2]);
			}
			if (defined(@entries[3])) { 
			  $featureset->signal(@entries[3]);
			}
			if (defined(@entries[4])) { 
			  $featureset->detection(@entries[4]);
			}
			if (defined(@entries[5])) { 
			  $featureset->detection_p-value(@entries[5]);
			}
			if (defined(@entries[6])) { 
			  $featureset->stat_common_pairs(@entries[6]);
			}
			if (defined(@entries[7])) { 
			  $featureset->signal_log_ratio(@entries[7]);
			}
			if (defined(@entries[8])) { 
			  $featureset->signal_log_ratio_low(@entries[8]);
			}
			if (defined(@entries[9])) { 
			  $featureset->signal_log_ratio_high(@entries[9]);
			}
			if (defined(@entries[10])) { 
			  $featureset->change(@entries[10]);
			}
			if (defined(@entries[11])) { 
			  $featureset->change_p-value(@entries[11]);
			}
			if (defined(@entries[12])) { 
			  $featureset->positive(@entries[12]);
			}
			if (defined(@entries[13])) { 
			  $featureset->negative(@entries[13]);
			}
			if (defined(@entries[14])) { 
			  $featureset->pairs(@entries[14]);
			}
			if (defined(@entries[15])) { 
			  $featureset->pairs_used(@entries[15]);
			}
			if (defined(@entries[16])) { 
			  $featureset->pairs_inavg(@entries[16]);
			}
			if (defined(@entries[17])) { 
			  $featureset->pos_fraction(@entries[17]);
			}
			if (defined(@entries[18])) { 
			  $featureset->log_avg(@entries[18]);
			}
			if (defined(@entries[19])) { 
			  $featureset->pos_neg(@entries[19]);
			}
			if (defined(@entries[20])) { 
			  $featureset->avg_diff(@entries[20]);
			}
			if (defined(@entries[21])) { 
			  $featureset->abs_call(@entries[21]);
			}
			if (defined(@entries[22])) { 
			  $featureset->inc(@entries[22]);
			}
			if (defined(@entries[23])) { 
			  $featureset->dec(@entries[23]);
			}
			if (defined(@entries[24])) { 
			  $featureset->inc_ratio(@entries[24]);
			}
			if (defined(@entries[25])) { 
			  $featureset->dec_ratio(@entries[25]);
			}
			if (defined(@entries[26])) { 
			  $featureset->pos_change(@entries[26]);
			}
			if (defined(@entries[27])) { 
			  $featureset->neg_change(@entries[27]);
			}
			if (defined(@entries[28])) { 
			  $featureset->inc_dec(@entries[28]);
			}
			if (defined(@entries[29])) { 
			  $featureset->dpos-dneg_ratio(@entries[29]);
			}
			if (defined(@entries[30])) { 
			  $featureset->log_avg_ratio_change(@entries[30]);
			}
			if (defined(@entries[31])) { 
			  $featureset->diff_call(@entries[31]);
			}
			if (defined(@entries[32])) { 
			  $featureset->avg_diff_change(@entries[32]);
			}
			if (defined(@entries[33])) { 
			  $featureset->b_a(@entries[33]);
			}
			if (defined(@entries[34])) { 
			  $featureset->fold_change(@entries[34]);
			}
			if (defined(@entries[35])) { 
			  $featureset->sort_score(@entries[35]);
			}
		  }
		}
	  }
	}
  }
}

1;
