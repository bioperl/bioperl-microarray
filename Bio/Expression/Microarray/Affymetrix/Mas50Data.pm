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

sub featuregroup {
  my($self,$arg) = @_;
  return $self->{featuregroup}->{$arg} if $self->{featureset}->{$arg};
  $self->{featuregroup}->{$arg} = Bio::Expression::FeatureGroup::FeatureGroupMas50->new()
	or $self->throw("Couldn't create a Bio::Expression::FeatureGroup::FeatureGroupMas50 $!");
  return $self->{featuregroup}->{$arg};
}

sub each_featuregroup {
  my $self = shift;
  my @return = ();
  foreach my $p (sort keys %{$self->{featuregroup}}) {
	push @return, $self->{featuregroup}->{$p};
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
		my $featuregroup = $self->array->featureset(@entries[1]);

		if (defined(@entries[2])) {
		  $featuregroup->probe_set_name(@entries[2]);
		}
		if (defined(@entries[3])) {
		  $featuregroup->stat_pairs(@entries[3]);
		}
		if (defined(@entries[4])) {
		  $featuregroup->stat_pairs_used(@entries[4]);
		}
		if (defined(@entries[5])) {
		  $featuregroup->signal(@entries[5]);
		}
		if (defined(@entries[6])) {
		  $featuregroup->detection(@entries[6]);
		}
		if (defined(@entries[7])) {
		  $featuregroup->detection_p_value(@entries[7]);
		}
		if (defined(@entries[8])) {
		  $featuregroup->stat_common_pairs(@entries[8]);
		}
		if (defined(@entries[9])) {
		  $featuregroup->signal_log_ratio(@entries[9]);
		}
		if (defined(@entries[10])) {
		  $featuregroup->signal_log_ratio_low(@entries[10]);
		}
		if (defined(@entries[11])) {
		  $featuregroup->signal_log_ratio_high(@entries[11]);
		}
		if (defined(@entries[12])) {
		  $featuregroup->change(@entries[12]);
		}
		if (defined(@entries[13])) {
		  $featuregroup->change_p_value(@entries[13]);
		}
		if (defined(@entries[14])) {
		  $featuregroup->positive(@entries[14]);
		}
		if (defined(@entries[15])) {
		  $featuregroup->negative(@entries[15]);
		}
		if (defined(@entries[16])) {
		  $featuregroup->pairs(@entries[16]);
		}
		if (defined(@entries[17])) {
		  $featuregroup->pairs_used(@entries[17]);
		}
		if (defined(@entries[18])) {
		  $featuregroup->pairs_inavg(@entries[18]);
		}
		if (defined(@entries[19])) {
		  $featuregroup->pos_fraction(@entries[19]);
		}
		if (defined(@entries[20])) {
		  $featuregroup->log_avg(@entries[20]);
		}
		if (defined(@entries[21])) {
		  $featuregroup->pos_neg(@entries[21]);
		}
		if (defined(@entries[22])) {
		  $featuregroup->avg_diff(@entries[22]);
		}
		if (defined(@entries[23])) {
		  $featuregroup->abs_call(@entries[23]);
		}
		if (defined(@entries[24])) {
		  $featuregroup->inc(@entries[24]);
		}
		if (defined(@entries[25])) {
		  $featuregroup->dec(@entries[25]);
		}
		if (defined(@entries[26])) {
		  $featuregroup->inc_ratio(@entries[26]);
		}
		if (defined(@entries[27])) {
		  $featuregroup->dec_ratio(@entries[27]);
		}
		if (defined(@entries[28])) {
		  $featuregroup->pos_change(@entries[28]);
		}
		if (defined(@entries[29])) {
		  $featuregroup->neg_change(@entries[29]);
		}
		if (defined(@entries[30])) {
		  $featuregroup->inc_dec(@entries[30]);
		}
		if (defined(@entries[31])) {
		  $featuregroup->dpos-dneg_ratio(@entries[31]);
		}
		if (defined(@entries[32])) {
		  $featuregroup->log_avg_ratio_change(@entries[32]);
		}
		if (defined(@entries[33])) {
		  $featuregroup->diff_call(@entries[33]);
		}
		if (defined(@entries[34])) {
		  $featuregroup->avg_diff_change(@entries[34]);
		}
		if (defined(@entries[35])) {
		  $featuregroup->b_a(@entries[35]);
		}
		if (defined(@entries[36])) {
		  $featuregroup->fold_change(@entries[36]);
		}
		if (defined(@entries[37])) {
		  $featuregroup->sort_score(@entries[37]);
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
			my $featuregroup = $self->array->featureset($self->array->name);

			if (defined(@entries[0])) { 
			  $featuregroup->probe_set_name(@entries[0]);
			}
			if (defined(@entries[1])) { 
			  $featuregroup->stat_pairs(@entries[1]);
			}
			if (defined(@entries[2])) { 
			  $featuregroup->stat_pairs_used(@entries[2]);
			}
			if (defined(@entries[3])) { 
			  $featuregroup->signal(@entries[3]);
			}
			if (defined(@entries[4])) { 
			  $featuregroup->detection(@entries[4]);
			}
			if (defined(@entries[5])) { 
			  $featuregroup->detection_p_value(@entries[5]);
			}
			if (defined(@entries[6])) { 
			  $featuregroup->stat_common_pairs(@entries[6]);
			}
			if (defined(@entries[7])) { 
			  $featuregroup->signal_log_ratio(@entries[7]);
			}
			if (defined(@entries[8])) { 
			  $featuregroup->signal_log_ratio_low(@entries[8]);
			}
			if (defined(@entries[9])) { 
			  $featuregroup->signal_log_ratio_high(@entries[9]);
			}
			if (defined(@entries[10])) { 
			  $featuregroup->change(@entries[10]);
			}
			if (defined(@entries[11])) { 
			  $featuregroup->change_p_value(@entries[11]);
			}
			if (defined(@entries[12])) { 
			  $featuregroup->positive(@entries[12]);
			}
			if (defined(@entries[13])) { 
			  $featuregroup->negative(@entries[13]);
			}
			if (defined(@entries[14])) { 
			  $featuregroup->pairs(@entries[14]);
			}
			if (defined(@entries[15])) { 
			  $featuregroup->pairs_used(@entries[15]);
			}
			if (defined(@entries[16])) { 
			  $featuregroup->pairs_inavg(@entries[16]);
			}
			if (defined(@entries[17])) { 
			  $featuregroup->pos_fraction(@entries[17]);
			}
			if (defined(@entries[18])) { 
			  $featuregroup->log_avg(@entries[18]);
			}
			if (defined(@entries[19])) { 
			  $featuregroup->pos_neg(@entries[19]);
			}
			if (defined(@entries[20])) { 
			  $featuregroup->avg_diff(@entries[20]);
			}
			if (defined(@entries[21])) { 
			  $featuregroup->abs_call(@entries[21]);
			}
			if (defined(@entries[22])) { 
			  $featuregroup->inc(@entries[22]);
			}
			if (defined(@entries[23])) { 
			  $featuregroup->dec(@entries[23]);
			}
			if (defined(@entries[24])) { 
			  $featuregroup->inc_ratio(@entries[24]);
			}
			if (defined(@entries[25])) { 
			  $featuregroup->dec_ratio(@entries[25]);
			}
			if (defined(@entries[26])) { 
			  $featuregroup->pos_change(@entries[26]);
			}
			if (defined(@entries[27])) { 
			  $featuregroup->neg_change(@entries[27]);
			}
			if (defined(@entries[28])) { 
			  $featuregroup->inc_dec(@entries[28]);
			}
			if (defined(@entries[29])) { 
			  $featuregroup->dpos-dneg_ratio(@entries[29]);
			}
			if (defined(@entries[30])) { 
			  $featuregroup->log_avg_ratio_change(@entries[30]);
			}
			if (defined(@entries[31])) { 
			  $featuregroup->diff_call(@entries[31]);
			}
			if (defined(@entries[32])) { 
			  $featuregroup->avg_diff_change(@entries[32]);
			}
			if (defined(@entries[33])) { 
			  $featuregroup->b_a(@entries[33]);
			}
			if (defined(@entries[34])) { 
			  $featuregroup->fold_change(@entries[34]);
			}
			if (defined(@entries[35])) { 
			  $featuregroup->sort_score(@entries[35]);
			}
		  }
		}
	  }
	}
  }
}

1;
