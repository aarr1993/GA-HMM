#!/usr/bin/env perl
use warnings;
use strict;
use Curses;

## Globals ##	
my $STARTLINE = 31;
my $curr_col = $STARTLINE; 
##

main();

sub main {
	initscr();

	my ($row, $col);
	getmaxyx($row, $col);

	my %hash = %{initHash()}; 
	printHeader($col, \%hash);

	## Start taking user input after end of header ##

	## mainloop ##
	
	while(1) {
		## Get user input ##

		my ($in) = userIn();

		## quit ## 

		if ($in eq "0") {
			last;
		}
	
		%hash = %{processUserIn($in, \%hash)};
		printHeader($col, \%hash);

		## Deal with wrapping lines ##

		if ($curr_col >= $row) {
			clearIn($row, $col, $STARTLINE);	
			$curr_col = $STARTLINE;
		}	
	}

	endwin();
}
	
sub initHash {
	my %hash;

	%{$hash{path}} = (
		wigpath => "/home/INIT",
		modelpath => "/home/INIT",
		countpath => "/home/INIT"
	);

	%{$hash{regions}} = (
		chrom => "INIT",
		run => [0,0],
		test => [0,0],
		model => [0,0]
	);

	%{$hash{params}} = (
		initialHMM => "INIT",
		templateHMM => "INIT",
		popSize => 0,
		gens => 0,
		outdir => ".",
		threads => 0,
		runfasta => "INIT"

	);
	
	%{$hash{opts}} = (
		norm => " ",
		customfa => " ",
		template => " ",
		initial => " ",
		cache => " ",
		ga => " "
	);	
	return (\%hash);
}

sub printHeader {
	my $col = $_[0];
	my %hash = %{$_[1]};
	my $str = "Path to wigs        : $hash{path}{wigpath}
Path to models      : $hash{path}{modelpath}
Path to count files : $hash{path}{countpath}

REGIONS
=======
chr   : [$hash{regions}{chrom}]
run   : [$hash{regions}{run}[0] : $hash{regions}{run}[1]]
test  : [$hash{regions}{test}[0] : $hash{regions}{test}[1]]
model : [$hash{regions}{model}[0] : $hash{regions}{model}[1]]

GA COMMAND
==========
./functionalGA.pl [run region fasta] [initial model] [template] [population size] [number of generations] [number of threads] [output dir]

./functionalGA.pl $hash{params}{runfasta} $hash{params}{initialHMM} $hash{params}{templateHMM} $hash{params}{popSize} $hash{params}{gens} $hash{params}{threads} $hash{params}{outdir}

#######################
# PREPROCESS COMMANDS #
#######################

0 Quit
1 Normalize wigfile      [$hash{opts}{norm}]  
2 Create .customfa files [$hash{opts}{customfa}]
3 Load template model    [$hash{opts}{template}]
4 Load initial model     [$hash{opts}{initial}]
5 Initialize cache       [$hash{opts}{cache}]
6 Prepare GA-HMM         [$hash{opts}{ga}]
7 Print batch file
";

	$str .= "#"x$col;
	addstring(0,0,$str);
	refresh();
}

sub userIn {
	addstring ($curr_col,0, "GA> ");
	refresh();
	my $input = getstring();
	$curr_col++;
	return $input;
}

sub clearIn {
	my $row = $_[0];
	my $col = $_[1];
#	my $start = $_[2];

	my $clear = " "x$col;

	for (my $i = $STARTLINE; $i <= $row; $i++) {
		addstring($i, 0, $clear);
	}

	refresh();
}

sub processUserIn {
	my $input = $_[0];
	my %hash = %{$_[1]};

	if ($input =~ /^[ch]/) {
		if ($input eq "h") {
			help();
		}
		elsif ($input =~ /^c/) {
			## run command ##
			my ($cmd) = $input =~ /c\s(.+)$/;
			my $cmd_out = system($cmd);
			## print output ##
		}
	} # letter commands: run bash command (c) and print help (h)
	elsif ($input =~ /[123456]/) {
		## call appropriate preprocessing function ##
		if ($input == 1) {
			%hash = %{normalize(\%hash)};
		} elsif ($input == 2) {
			%hash = %{customfa(\%hash)};
		} elsif ($input == 3) {
			%hash = %{templateModel(\%hash)};
		} elsif ($input == 4) {
			%hash = %{initialModel(\%hash)};
		} elsif ($input == 5) {
			%hash = %{initCache(\%hash)};
		} elsif ($input == 6) {
			%hash = %{runGA(\%hash)};
		}
		elsif ($input == 7) {
			printbatch(\%hash);
		}
	} # input is a preprocessing command
	else {
		error(0); #invalid input error
	}
	return(\%hash);
}

sub text {
	## prints text to console and cleans up stuff ##
	my $str = $_[0];
	my ($y, $x);
	getmaxyx($y, $x);
	my $usedrows = int( length($str) / $x) + 1;

	if ($curr_col + $usedrows >= $y) {
		clearIn($y, $x, $STARTLINE);
		$curr_col = $STARTLINE;
	}
	
	addstring($curr_col, 0, $str);
	$curr_col += $usedrows;
}

sub error {
	## deal with any errors ##
	my $err = $_[0];
	
	if ($err == 0) {
		## print invalid input error ##
		text("ERROR 0: Invalid input");
	}
	elsif ($err == 1) {
		text("ERROR 1: Insufficient prereqs to execute function. Function is [] prereqs are []");
	}
}

sub help {
	## print help guide and usage ##
}

sub normalize {
	my %hash = %{$_[0]};
	my $in;

	text("Enter chr wigfile (/path/to/wig)");	
	my $wig = userIn();
	text("Is $wig normalized? (y/n)");
	$in = userIn();
	
	if ($in eq "y") {
		text("Does the default [../wigs/chr_wig_shifted.wig] have the normalized wig? (y/n)");
		$in = userIn();
		if ($in eq 'y') {
			$hash{path}{wigpath} = "../wigs/chr_wig_shifted.wig";
			$hash{opts}{norm} = "X";
		}
		elsif ($in eq 'n') {
			$hash{path}{wigpath} = $wig;
			$hash{opts}{norm} = "X";
		}
		else {
			text("Bad input! Try again.");
		}
	}
	elsif ($in eq 'n') {
		## normalize wig ##
		text("Normalizing");
		system("../subscripts/normalize.pl $wig ../wigs/chr_wig_shifted.wig");
		text("Normalized wig is in ../wigs/chr_wig_shifted.wig");
		$hash{path}{wigpath} = "../wigs/chr_wig_shifted.wig";
		$hash{opts}{norm} = "X";
	}
	else {
		text("Bad input! Try again.");
	}
		
	return (\%hash);
}	

sub customfa {
	my %hash = %{$_[0]};
	my $in;
	my ($start, $end);

	text("Enter chromosome.");
	$hash{regions}{chrom} = userIn();

	text("Enter space-seperated model region [ex. 100 200].");
	$in = userIn();
	($start, $end) = split(/\s/, $in);
	$hash{regions}{model}[0] = $start;
	$hash{regions}{model}[1] = $end;

	text("Enter space-seperated run region [ex. 100 200].");
	$in = userIn();
	($start, $end) = split(/\s/, $in);
	$hash{regions}{run}[0] = $start;
	$hash{regions}{run}[1] = $end;

	text("Enter space-seperated test region (within run region) [ex. 100 200].");
	$in = userIn();
	($start, $end) = split(/\s/, $in);
	$hash{regions}{test}[0] = $start;
	$hash{regions}{test}[1] = $end;
	
	text(".customfa files and [run, model, test] wigs already created? (y/n)");
	$in = userIn();

	if ($in eq 'y') {
		text("Assuming default filenames: $hash{path}{wigpath} ../fasta/chr.customfa ../{wigs, fasta}/run_region.{wig, customfa} and ../{wigs, fasta}/test_region.{wig, customfa}");
		
		$hash{params}{runfasta} = "../fasta/run.customfa";
		$hash{opts}{customfa} = "X";
	} 
	elsif ($in eq 'n') {
		text("Making small wigs and customfa files");
		# test
		system("../subscripts/extract_wig.pl $hash{path}{wigpath} $hash{regions}{chrom} $hash{regions}{test}[0] $hash{regions}{test}[1] ../wigs/test_region.wig");
		system("../subscripts/wig2fa.pl -i ../wigs/test_region.wig -o ../fasta/test_region.customfa");
		# run
		system("../subscripts/extract_wig.pl $hash{path}{wigpath} $hash{regions}{chrom} $hash{regions}{run}[0] $hash{regions}{run}[1] ../wigs/run_region.wig");
		system("../subscripts/wig2fa.pl -i ../wigs/run_region.wig -o ../fasta/run_region.customfa");
		# chrom
		system("../subscripts/wig2fa.pl -i $hash{path}{wigpath} -o ../fasta/chr.customfa");
		$hash{opts}{customfa} = "X";
		$hash{params}{runfasta} = "../fasta/run.customfa";
				
	} else {
		text("Bad input! Try again.");
	}		
	return(\%hash);
}

sub templateModel {
	my %hash = %{$_[0]};
	my $in;
	
	text("Enter the path to directory of model region beds.");
	my $path = userIn();
	
	$hash{path}{modelpath} = $path;
	
	text("Template HMM created? (y/n)");
	$in = userIn();
	
	if ($in eq "y") {
		text("Enter the name of template.");
		my $template = userIn();
		$hash{params}{templateHMM} = $template;
		$hash{opts}{template} = "X";
	}
	elsif ($in eq "n") {
		%hash = %{create_template(\%hash)};
		$hash{opts}{template} = "X";
	}
	else {
		text("Bad input! Try again.");
	}
	return(\%hash);
}	

sub initialModel {
	my %hash = %{$_[0]};

# needs template HMM to be created first
	text("Enter filename for initial model");
	my $name = userIn();
	open (DEBUG, ">", "err.log");

	if ($hash{path}{modelpath} !~ /INIT/) {
		print DEBUG "[$hash{path}{modelpath}]\n";		
		open (TEMP, "<", $hash{params}{templateHMM}) or (text("Make template HMM first!") && return(\%hash));
		

		my $index = 0;
		my @sections;

		while (<TEMP>) {
			my $line = $_;
			

			if ($line =~ /\#{45}/) {
				$index++;
			}

			$sections[$index] .= $line;			

		}
		close TEMP;
		print DEBUG "Sections are $index\n";
		my ($order) = $sections[2] =~ /ORDER:\t(\d)/;
		my ($emmisions) = $sections[0] =~ /SCORE:\s(.+)/;	
		print DEBUG "emmissions are $emmisions\n";
		close DEBUG;
		my @files = <$hash{path}{modelpath}*>;
		
		$index = 2;
		open (MODEL, ">", "../models/$name");
		print MODEL $sections[0];
		print MODEL $sections[1];


		for (my $i = 0; $i < @files; $i++) {
			my $dir;
			($dir, $files[$i]) = $files[$i] =~ /^(.+)\/([^\/]+)$/;

			system("fastaFromBed -fi ../fasta/chr.customfa -bed $dir\/$files[$i] -fo ../fasta/$files[$i].fa");
			system("../subscripts/HMM_Counter.pl -i ../fasta/$files[$i].fa -r $order -w $emmisions -o ../emm/$files[$i].count");
			$hash{path}{countpath} = "../emm/";
			
			my $emm;

			open (BED, "<", "$hash{path}{countpath}$files[$i].count");
			while(<BED>) {
				my $line = $_;
				$emm .= $line;
			}
			close BED;

			print MODEL $emm;
		}

		print MODEL $sections[@sections - 1];
		close MODEL;
		
		$hash{params}{initialHMM} = "../models/$name";
		$hash{opts}{initial} = "X";
	}
	else {
		text("Create template HMM first");
	}

	return(\%hash);
}

sub create_template {
	my %hash = %{$_[0]};

	text("Creates a 7 state model, with 3 NOISY and 3 SPARSE states, and 2 peak states per.");

	text("Enter model name, no spaces.");
	my $name = userIn();
	my $date = localtime();
	my @files = <$hash{path}{modelpath}*>;
	
	text("Using filenames as peak state names. Format should be [noisy|sparse|genomic]_[inter|broad|med|sharp]peak.bed Print example? (y/n)");
	my $in = userIn();

	if ($in eq 'y') {
		text("Ex. noisy_broadpeak.bed noisy_sparsepeak.bed noisy_interpeak.bed sparse_medpeak.bed sparse_sharppeak.bed sparse_interpeak.bed genomic_interpeak.bed");
	}

	text("Processing files ...");
	
	if (@files != 7) {
		text("Error: Number of files != 7!");
	}	

	my @states;

	for (my $i = 0; $i < @files; $i++) {
		my ($state_name) = $files[$i] =~ /([^\/]+)\.bed$/;
		($state_name) =~ tr/[a-z]/[A-Z]/;

		$states[$i]{name} = $state_name;
	
		open (IN, "<", $files[$i]) or text("Cannot open $files[$i]");
		
		my $avg = 0;
		my $count = 0;

		while(<IN>) {
			my $line = $_;
			chomp $line;
			my ($chr, $start, $end) = split(/\t/, $line);
			$avg += ($end - $start);
			$count++;
		}	
		$avg = (1/($avg/$count));
		
		## setting transition probabilities ##
		if ($state_name !~ /INTER/) {
			my $num = (1 - $avg)/3;
			if ($state_name =~ /SPARSE/) {
				$states[$i]{trans}{NOISY_BROADPEAK} = 0;
				$states[$i]{trans}{NOISY_MEDPEAK} = 0;
				$states[$i]{trans}{NOISY_INTERPEAK} = 0;
				$states[$i]{trans}{GENOMIC_INTERPEAK} = $num;
				$states[$i]{trans}{SPARSE_INTERPEAK} = $num*2;
				if ($state_name =~ /SHARP/) {
					$states[$i]{trans}{SPARSE_SHARPPEAK} = $avg;
					$states[$i]{trans}{SPARSE_MEDPEAK} = 0;
				} 
				else {
					$states[$i]{trans}{SPARSE_SHARPPEAK} = 0;
					$states[$i]{trans}{SPARSE_MEDPEAK} = $avg;
				}
			}
			else {
				$states[$i]{trans}{SPARSE_SHARPPEAK} = 0;
				$states[$i]{trans}{SPARSE_MEDPEAK} = 0;
				$states[$i]{trans}{SPARSE_INTERPEAK} = 0;
				$states[$i]{trans}{GENOMIC_INTERPEAK} = $num;
				$states[$i]{trans}{NOISY_INTERPEAK} = $num*2;
				if ($state_name =~ /BROAD/) {
					$states[$i]{trans}{NOISY_BROADPEAK} = $avg;
					$states[$i]{trans}{NOISY_MEDPEAK} = 0;
				} 
				else {
					$states[$i]{trans}{NOISY_BROADPEAK} = 0;
					$states[$i]{trans}{NOISY_MEDPEAK} = $avg;
				}

			}
		}
		else {
			if ($state_name =~ /GENOMIC/) {
				my $num = (1 - $avg) / 4;
				$states[$i]{trans}{NOISY_BROADPEAK} = $num;
				$states[$i]{trans}{NOISY_MEDPEAK} = $num;
				$states[$i]{trans}{SPARSE_SHARPPEAK} = $num;
				$states[$i]{trans}{SPARSE_MEDPEAK} = $num;
				$states[$i]{trans}{GENOMIC_INTERPEAK} = $avg;
				$states[$i]{trans}{NOISY_INTERPEAK} = 0;
				$states[$i]{trans}{SPARSE_INTERPEAK} = 0;
			}
			else {
				my $num = (1 - $avg) / 2;

				if ($state_name =~ /NOISY/) {
					$states[$i]{trans}{NOISY_BROADPEAK} = $num;
					$states[$i]{trans}{NOISY_MEDPEAK} = $num;
					$states[$i]{trans}{NOISY_INTERPEAK} = $avg;
					$states[$i]{trans}{SPARSE_SHARPPEAK} = 0;
					$states[$i]{trans}{SPARSE_MEDPEAK} = 0;
					$states[$i]{trans}{SPARSE_INTERPEAK} = 0;
					$states[$i]{trans}{GENOMIC_INTERPEAK} = 0;
				}
				else {
					$states[$i]{trans}{SPARSE_SHARPPEAK} = $num;
					$states[$i]{trans}{SPARSE_MEDPEAK} = $num;
					$states[$i]{trans}{SPARSE_INTERPEAK} = $avg;
					$states[$i]{trans}{NOISY_BROADPEAK} = 0;
					$states[$i]{trans}{NOISY_MEDPEAK} = 0;
					$states[$i]{trans}{NOISY_INTERPEAK} = 0;
					$states[$i]{trans}{GENOMIC_INTERPEAK} = 0;
				}
			}
		}		

		close IN;
	}
	open (OUT, ">", "../models/$name.template");
	
	print OUT "#STOCHHMM MODEL FILE\nMODEL INFORMATION\n" . "="x54 . "\n";
	print OUT "MODEL_NAME:\t$name\nMODEL_DESCRIPTION:\t\nMODEL_CREATION_DATE:\t$date\n\nTRACK SYMBOL DEFINITIONS\n" . "="x54 . "\nSCORE: N,L,O,M,H,S\n\nSTATE DEFINITIONS\n" . "\#"x45 . "\n";

	print OUT "STATE:\n\tNAME:\tINIT\nTRANSITION:\tSTANDARD:\tP(X)\n";

	foreach my $tr (keys %{$states[0]{trans}}) {
		print OUT "\t$tr\t" . (1/7) . "\n"; 
	}
		
	print OUT "\#"x45 . "\n";

	my $letter = "A";

	for (my $j = 0; $j < @states; $j++) {
		print OUT "STATE:\n\tNAME:\t$states[$j]{name}\n\tPATH_LABEL:\t$letter\n\tGFF_DESC:$states[$j]{name}\nTRANSITION:\tSTANDARD:\tP(X)\n";

		foreach my $transition (keys %{$states[$j]{trans}}) {
			print OUT "\t$transition\t$states[$j]{trans}{$transition}\n";
		}

		print OUT "\tEND:\t1\nEMISSION:\tSCORE\tCOUNTS\n\tORDER:\t3\n\@N\tL\tO\tM\tH\tS\n";
		print OUT "\#"x45 . "\n";
	
		$letter++;
	}

	print OUT "//END";
	close OUT;
	
	$hash{params}{templateHMM} = "../models/$name.template";
	return (\%hash);
}

sub initCache {
	my %hash = %{$_[0]};

	text("Initializing cache for eval script");
	my ($blocks) = `"../subscripts/make_cache.pl ../cache ../wigs/test_region.wig"`; # backticks necessary to capture command output
	text("blocks is $blocks");
	open (TMP, ">", "../tmp/sig_unsig_blocks.txt");
	print TMP "test_region.wig\t../cache\t$blocks\n";
	close TMP;

	$hash{opts}{cache} = "X";

	return (\%hash);
}

sub runGA {
	my %hash = %{$_[0]};

	text("Enter population size");
	my $pop = userIn();
	text("Enter number of generations");
	my $gens = userIn();	
	text("Enter number of threads");
	my $threads = userIn();
	text("Enter output directory");
	my $out = userIn();

	$hash{params}{popSize} = $pop;	
	$hash{params}{gens} = $gens;	
	$hash{params}{threads} = $threads;	
	$hash{params}{outdir} = $out;

	$hash{opts}{ga} = "X";

	return (\%hash);
}

sub printbatch {
	my %hash = %{$_[0]};
	
	text("Enter name of batch file");
	my $name = userIn();
	text("Enter name of job");
	my $job = userIn();
	text("Enter email for updates");
	my $email = userIn();

	my $batch = "#!/bin/bash -l
# NOTE the -l flag!

# If you need any help, please email help\@cse.ucdavis.edu

# Name of the job - You'll probably want to customize this.
#SBATCH -J $job

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o $job.output
#SBATCH -e $job.output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$email

# no -n here, the user is expected to provide that on the command line.

# The useful part of your job goes below

# run one thread for each one the user asks the queue for
# hostname is just for debugging
hostname
export OMP_NUM_THREADS=\$SLURM_NTASKS

./functionalGA.pl ../fasta/run_region.customfa $hash{params}{initialHMM} $hash{params}{templateHMM} $hash{params}{popSize} $hash{params}{gens} $hash{params}{threads} $hash{params}{out}
";

	open (BATCH, "<", $name);
	print BATCH $batch;
	close BATCH;
	
	text("Printed batch file to $batch. Please run it to run GA. Threads is $hash{params}{threads}");
}
