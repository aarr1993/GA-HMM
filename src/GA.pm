package GA;

use strict;
use warnings;

our (%wigfiles, %regions, $chromosome, $template_hmm, $initial_hmm, $testvar);

## Initialization ##

%wigfiles = (
  chr_wig => "",
  test_wig => "",
  run_wig => "",
); 

%regions = (
  model_region => [0,0], # is this needed?
  test_region => [0,0],
  run_region => [0,0],
);

##

sub test {
  print "GA.pm test\n";
  print "uninitialized test var is $testvar\n";
  $testvar = 20;
  print "testvar initialized to 20 is $testvar\n";
} 

sub normalize {
  print "Enter entire chr wig filename: ";
  my $chr_wig = <>;
  chomp $chr_wig;
  
  $wigfiles{chr_wig} = $chr_wig;

  print "Already normalized? (y/n): ";
  my $check = <>;
  chomp $check;

  return if ($check =~ /y/i) ;
  # trusting user on file format

  # linear shift based on > 3rd quartile median.
  `../subscripts/normalize.pl $wigfiles{chr_wig} ../wigs/chr_wig_shifted.wig`;
  $wigfiles{chr_wig} = "../wigs/chr_wig_shifted.wig";
}

sub customfa {
  print ".customfa files and small wigs already created? (y/n): ";
  my $in = <>;
  chomp $in;

  if ($in =~ /y/i) {
    print "Using default filenames: ../wigs/run_region.wig and ../wigs/test_region.wig\n";
    $wigfiles{run_wig} = '../wigs/run_region.wig';
    $wigfiles{test_wig} = '../wigs/test_region.wig';
    return;
  } 

  print "Initialize regions from file? (y/n): ";
  $in = <>;
  chomp $in;

  if ($in =~ /y/i) {
    print "File should be in format:
chr8
MODEL 3000000 14000000
TEST 15000000 42000000
RUN 15000000 120000000

Enter filename: ";
    my $file = <>;
    chomp $file;
    open (IN, "<", $file) or die "Could not open file\n";
    while(<IN>) {
      my $line = $_;
      chomp $line;
      if ($line =~ /^chr/) {
        $chromosome = $line =~ /^(.+)$/; 
      }
      elsif ($line =~ /^MODEL/) {
        my @tmp = $line =~ /^MODEL\s+(\d+)\s+(\d+)/;
        $regions{model_region} = \@tmp;
      }
      elsif ($line =~ /TEST/) {
        my @tmp = $line =~ /^TEST\s+(\d+)\s+(\d+)/;
        $regions{test_region} = \@tmp;
      }
      elsif ($line =~ /RUN/) {
        my @tmp = $line =~ /^RUN\s+(\d+)\s+(\d+)/;
        $regions{run_region} = \@tmp;
      }
    }
    close IN;
    print "Regions read from file $file\n";
  }  
  else {
    print "Enter chromosome used for this processs (ex. chr8): ";
    my $in_chr = <>;
    chomp $in_chr;
    $chromosome = $in_chr;

    print "Enter coordinates for the model region, space seperated, no commas: ";
    my $input = <>;
    chomp $input;
    my @model = split(/\s/, $input);
    $regions{model_region} = \@model; # FIXME correct syntax???

    print "Enter coordinates for test region: ";
    $input = <>;
    chomp $input;
    my @test = split(/\s/, $input);
    $regions{test_region} = \@test;

    print "Enter coordinates for run region: ";
    $input = <>;
    chomp $input;
    my @run = split(/\s/, $input);
    $regions{run_region} = \@run;
  }

  ## User input done ##
  
  print "Outputting small region wigs and making .customfa files\n";

  # Test region 
   `../subscripts/extract_wig.pl $wigfiles{chr_wig} $chromosome $regions{test_region}[0] $regions{test_region}[1] ../wigs/test_region.wig` ;
  $wigfiles{test_region} = '../wigs/test_region.wig';
  `../subscripts/wig2fa.pl -i ../wigs/test_region.wig -o ../fasta/test_region.customfa`;

  # Run region
   `../subscripts/extract_wig.pl $wigfiles{chr_wig} $chromosome $regions{run_region}[0] $regions{run_region}[1] ../wigs/run_region.wig` ;
  $wigfiles{run_region} = '../wigs/run_region.wig';
  `../subscripts/wig2fa.pl -i ../wigs/run_region.wig -o ../fasta/run_region.customfa` ;

  # Entire chromosome
  `../subscripts/wig2fa.pl -i $wigfiles{chr_wig} -o ../fasta/chr.customfa` ;

  # TODO
  # die if error
}

sub set_template {
  my $check = "n";

  do {
    print "Template HMM already created? (y/n): ";
    $check = <>;
    chomp $check;

    if ($check =~ /y/i) {
      print "Enter filename of template HMM: ";
      my $file = <>;
      chomp $file;
      $template_hmm = $file;
    }
    else {
      print "Please open new terminal and run generate_HMM script. This option will loop until 'y' is input.\n";
    }   
  } while ($check !~ /y/i);
}

sub emissions {
  my @beds;

  my $check = "n";
  my $dir = '';
  do {
    print "Bedfiles for initial hmm already created? (y/n): ";
    $check = <>;
    chomp $check;
  
    if ($check =~ /y/i) {
      print "Enter directory path to beds: ";
      $dir = <>;
      chomp $dir;
      print "Enter space-seperated filenames: ";
      my $input = <>;
      chomp $input;
      @beds = split(/\s/, $input);
    } 
    else {
      print "Within the region [" . $regions{model_region}[0] . "," . $regions{model_region}[1] . "] find small regions that correspond to each state in the template HMM file. Make a seperate bedgraph file (chr\tstart\tend) for each state. Select 3-4 small regions per state. This option will loop until 'y' is input.\nIt is recommended to use the UCSC Genome Browser or other method to view the wig file and manually select regions that correspond to the states.\n";
    }
  } while ($check !~ /y/i);

  print "Creating .count files for initial HMM emissions.\n";

  print "Enter order: ";
  my $order = <>;
  chomp $order;
  print "Enter emissions comma-seperated EXACTLY as in template/initial HMM model file: ";
  my $emm = <>;
  chomp $emm;

  for (my $i = 0; $i < @beds; $i++) {
    `fastaFromBed -fi ../fasta/chr.customfa -bed $dir/$beds[$i] -fo ../fasta/$beds[$i].fa`;
    `../subscripts/HMM_Counter.pl -i ../fasta/$beds[$i].fa -r $order -w $emm -o ../emm/$beds[$i].count`;
  }

  print ".count files created. They are in the format <bedfile>.count in the emm/ directory. Please append these count files beneath the correct state in the 'initial HMM' file.\n";
  print "Once this is done, input filename of initial HMM file: ";
  my $file = <>;
  chomp $file;
  $initial_hmm = $file; 
  # FIXME automate this?
}

sub set_model {
  if ($initial_hmm) {
    print "initial HMM file is already set to $initial_hmm. Go on to next step\n";
  }
  else {
    print "Input filename of initial HMM file: ";
    my $file = <>;
    chomp $file;
    $initial_hmm = $file; 
  }
}

sub initialize {
  print "Initializing cache for evaluation script with $wigfiles{test_wig}.\n";
  
  my $blocks = `../subscripts/make_cache.pl ../cache ../wigs/test_region.wig`;
  open (TMP, ">", "../tmp/sig_unsig_blocks.txt");
  print TMP "test_region.wig\t../cache\t$blocks\n";
  close TMP;
  # eval script will read from this file for its configuration at each call
  # may slow it down
}

sub params {
  print "Setting parameters to run GA\n";
  print "run region customfa is ../fasta/run_region.customfa\n"; 
  print "initial hmm is $initial_hmm\n";
  print "template hmm is $template_hmm\n";

  print "Enter population size: ";
  my $popsize = <>;
  chomp $popsize;

  print "Enter number of generations: ";
  my $gens = <>;
  chomp $gens;

  print "Enter number of threads: ";
  my $threads = <>;
  chomp $threads;

  print "Enter output directory/ : ";
  my $outdir = <>;
  chomp $outdir;

  my $cmd = "./functionalGA.pl ../fasta/run_region.customfa $initial_hmm $template_hmm $popsize $gens $threads $outdir";

  print "Command is\n" . $cmd . "\n";
  print "Run GA? (y/n): ";
  my $option = <>;
  chomp $option;

  # TODO
  # do some sort of checks / file sanity stuff here  

  if ($option =~ /y/i) {
    `$cmd`;
  }
}

### User validation functions ###

# validate_file_type($filename, ".extension");
sub validate_file_type {
  my $filename = $_[0];
  my $ext = $_[1];
  
  if (! -f $filename) {
    die "File does not exist at $filename!\n";
  }

  if ($filename !~ /$ext$/) {
    die "Incorrect file type! Was expecting $ext, got $filename\n";
  } 
}

# validate_wig($filename);
sub validate_wig {
  my $filename = $_[0];

  validate_file_type($filename, ".wig");
  open (IN, "<", $filename) or die "Could not open file $filename\n"; # also checks that it exists

  my $prev_span = -1;
  my $change = -1;

  while(<IN>) {
    my $line = $_;
    chomp $line;
    
    if ($line =~ /Step/) {
      my ($span) = $line =~ /Step chrom=.+ span=(\d+)/;
      if ($span != $prev_span) {
        $change++;
      }  
    }     
  }
  close IN;

  if ($change) {
    die "span in $filename changed $change times! Span MUST be constant in wig file!\n";
  }
}

# validate_regions();
sub validate_regions {
  foreach my $field (keys %regions) {
    if ($regions{$field}[0] == 0 && $regions{$field}[1] == 0) {
      die "region $field coordinates are not set!\n";
    }
    elsif ($regions{$field}[0] == $regions{$field}[1]) {
      die "region $field start and end are the same!\n";
    }
    elsif ($regions{$field}[1] < $regions{$field}[0]) {
      die "region $field end < start!";
    }
  }
  if ($regions{model_region}[1] > $regions{run_region}[0] && $regions{model_region}[0] < $regions{run_region}[1]) {
    die "Model region is within run region!\n";
  }
  if (!($regions{test_region}[1] >= $regions{run_region}[0] && $regions{test_region}[0] <= $regions{run_region}[1])) {
    die "Test region must be within run region!\n";
  }
}

# validate_fasta($filename);
sub validate_fasta {
  my $filename = $_[0];

}

# validate_template_hmm($filename);
sub validate_template_hmm {
  my $filename = $_[0];

}
