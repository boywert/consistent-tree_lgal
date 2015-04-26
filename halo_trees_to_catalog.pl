#!/usr/bin/perl -w
use lib qw(src);
use ReadConfig;
use MultiThread;
use Universe::Time;
use MassConversions;
use IO::File;
use Scalar::Util;

$MultiThread::threadlimit = 20;
load_config();

opendir DIR, $TREE_OUTBASE;
my @trees = grep { /^tree_.*\.dat$/ } readdir DIR;
closedir DIR;

open INPUT, "<", "$TREE_OUTBASE/$trees[0]";
our $firstline = <INPUT>;
chomp($firstline);
my @elems = split(" ", $firstline);
push @elems, qw/Macc Mpeak Vacc Vpeak Halfmass_Scale Acc_Rate_Inst Acc_Rate_100Myr Acc_Rate_1*Tdyn Acc_Rate_2*Tdyn Acc_Rate_Mpeak Mpeak_Scale Acc_Scale First_Acc_Scale First_Acc_Mvir First_Acc_Vmax Vmax@Mpeak/;

for (0..$#elems) {
    $elems[$_] =~ s/\(\d+\)$//;
    $elems[$_].="($_)";
}
$firstline = "@elems\n";
while (<INPUT>) {
    last unless (/^#/);
    $firstline .= $_;
}
close INPUT;
$firstline .= "#Macc,Vacc: Mass and Vmax at accretion.\n";
$firstline .= "#Mpeak,Vpeak: Peak mass and Vmax over mass accretion history.\n";
$firstline .= "#Halfmass_Scale: Scale factor at which the MMP reaches 0.5*Mpeak.\n";
$firstline .= "#Acc_Rate_*: Halo mass accretion rates in Msun/h/yr.\n";
$firstline .= "#            Inst: instantaneous; 100Myr: averaged over past 100Myr,\n";
$firstline .= "#            X*Tdyn: averaged over past X*virial dynamical time.\n";
$firstline .= "#            Mpeak: Growth Rate of Mpeak, averaged from current z to z+0.5\n";
$firstline .= "#Mpeak_Scale: Scale at which Mpeak was reached.\n";
$firstline .= "#Acc_Scale: Scale at which satellites were (last) accreted.\n";
$firstline .= "#First_Acc_Scale: Scale at which current and former satellites first passed through a larger halo.\n";
$firstline .= "#First_Acc_(Mvir|Vmax): Mvir and Vmax at First_Acc_Scale.\n";
$firstline .= "#Vmax\@Mpeak: Halo Vmax at the scale at which Mpeak was reached.\n";


opendir DIR, $HLIST_OUTBASE;
for (readdir DIR) {
    next unless /^hlist_[0-9.]+\.list$/;
    unlink "$HLIST_OUTBASE/$_";
}
closedir DIR;

foreach (@trees) {
    next unless (MultiThread::ForkChild());
    $main::filename = $_;
    convert_to_catalog("$TREE_OUTBASE/$_");
    exit;
}
MultiThread::WaitForChildren();

opendir DIR, $HLIST_OUTBASE;
for (readdir DIR) {
    next unless /^hlist_([0-9.]+)\.list\.tree/;
    my $scale = $1;
    push @{$files{$scale}}, $_;
}
closedir DIR;

foreach my $scale (keys %files) {
    next unless (MultiThread::ForkChild());
    my @to_combine = @{$files{$scale}};
    chdir $HLIST_OUTBASE;
    open OUTPUT, ">", "hlist_$scale.list" or 
	die "Couldn't open output file hlist_$scale.list! [$!]\n";;
    print OUTPUT $firstline;
    my $buffer;
    for (sort @to_combine) {
	open INPUT, "<", $_ or die "Couldn't open input file $_! [$!]\n";
	while (($n = read(INPUT, $buffer, 1000000))) {
	    print OUTPUT $buffer;
	}
	die "Failed to read from $_ [$!]!\n" unless (defined $n);
    }
    close OUTPUT;
    unlink @to_combine;
    exit;
}
MultiThread::WaitForChildren();


our (%times, %dyn_times);

sub convert_to_catalog 
{
    @ARGV = shift;
    while (<>) {
	next unless /\#\s*tree/i; #Skip header
	last;
    }

    while (<>) {
	if (/^\#\s*tree/) {
	    process_tree();
	    %halos = ();
	    @halos = ();
	    next;
	}
	next unless (/^\s*\d+\.?\d*\s+\d+/);
	my $h = new TreeHalo;
	$h->scanline($_);
	push @halos, $h;
	$halos{$h->{scale}}{$h->{id}} = $h;
	if (!$times{$h->{scale}}) {
	    $times{$h->{scale}} = scale_to_years($h->{scale});
	    $dyn_times{$h->{scale}} =
		MassConversions::Dyn_time(1.0/$h->{scale} - 1, $Om, $h0);
	    $zphalf_times{$h->{scale}} = scale_to_years($h->{scale}) - scale_to_years(1.0/(0.5+1.0/$h->{scale}));
	}
    }
    process_tree();
    $_->close() foreach (values %TreeHalo::tree_outputs);
}

sub max {
    return (($_[0] > $_[1]) ? $_[0] : $_[1]);
}

sub process_tree {
    my $h;
    for $h (@halos) {
	$h->{num_prog} = 0;
    }
    for $h (@halos) {
	my $d = $halos{$h->{desc_scale}}{$h->{descid}};
	next unless defined($d);
	$d->{num_prog}++;
	my $p = $d->{prog};
	if (defined($p) and $p->{mvir} > $h->{mvir}) {
	    next;
	}
	$d->{prog} = $h;
    }
    for $h (@halos) {
	calc_mass_vmax_acc($h);
	calc_halfmass($h);
	calc_accretion_rates($h);
	$h->print();
    }
}

sub calc_mass_vmax_acc {
    no warnings 'recursion';
    my $h = shift;
    return if ($h->{seen});
    $h->{seen} = 1;
    if ($h->{prog}) {
	calc_mass_vmax_acc($h->{prog});
	if ($h->{upid}>-1) { # If we are a subhalo
	    $h->{vacc} = $h->{prog}{vacc};
	    $h->{macc} = $h->{prog}{macc};
	    $h->{acc_scale} = $h->{prog}{acc_scale};
	    $h->{first_acc} = $h->{prog}{first_acc};
	} else {
	    $h->{vacc} = $h->{vmax};
	    $h->{macc} = $h->{orig_mvir};
	    $h->{acc_scale} = $h->{scale};
	    $h->{first_acc} = $h;
	    my $hpf = $h->{prog}{first_acc};
	    $h->{first_acc} = $hpf if ($hpf and ($hpf->{id} != $h->{prog}{id}) and $hpf->{mpeak}*2.0 > $h->{orig_mvir});
	}
	$h->{vpeak} = max($h->{vmax}, $h->{prog}{vpeak});

	$h->{mpeak_scale} = $h->{prog}{mpeak_scale};
	$h->{vmpeak} = $h->{prog}{vmpeak};
	$h->{mpeak} = $h->{prog}{mpeak};
	if ($h->{orig_mvir} > $h->{prog}{mpeak}) {
	    $h->{mpeak} = $h->{orig_mvir};
	    $h->{mpeak_scale} = $h->{scale};
	    $h->{vmpeak} = $h->{vmax};
	}

	#Vpeak / Mpeak *before* accretion
	#$h->{vpeak} = max($h->{vacc}, $h->{prog}{vpeak});
	#$h->{mpeak} = max($h->{macc}, $h->{prog}{mpeak});
    } else {
	$h->{vpeak} = $h->{vmax};
	$h->{vmpeak} = $h->{vmax};
	$h->{vacc} = $h->{vmax};
	$h->{mpeak} = $h->{orig_mvir};
	$h->{macc} = $h->{orig_mvir};
	$h->{acc_scale} = $h->{scale};
	$h->{mpeak_scale} = $h->{scale};
	$h->{first_acc} = $h;
    }
    Scalar::Util::weaken($h->{first_acc});
}

sub calc_halfmass {
    no warnings 'recursion';
    my $h = shift;
    return if ($h->{seen}>1);
    $h->{seen} = 2;
    if ($h->{prog}) {
	my $hm = calc_halfmass($h->{prog});
	my $hm2 = $hm;
	while (defined($hm2 = $halos{$hm2->{desc_scale}}{$hm2->{descid}})) {
	    last unless ($hm2->{orig_mvir} < 0.5*$h->{mpeak});
	    $hm = $hm2 if ($hm2->{orig_mvir} > $hm->{orig_mvir});
	}
	$h->{halfmass} = $hm->{scale};
	return($hm);
    }
    else {
	$h->{halfmass} = $h->{scale};
	return $h;
    }
}

sub find_mass_years_ago {
    my $h = shift;
    my $guess = shift;
    my $t = shift;
    my $tnow = $times{$h->{scale}};
    my $tthen = $times{$guess->{scale}};
    while (defined ($hm2 = $halos{$guess->{desc_scale}}{$guess->{descid}}) and
	   ($tnow-$tthen) > $t) {
	$guess = $hm2;
	$tthen = $times{$hm2->{scale}};
    }

    while ($guess->{prog} and ($tnow - $tthen < $t)) {
	$guess = $guess->{prog};
	$tthen = $times{$guess->{scale}};
    }

    my $guess2 = $halos{$guess->{desc_scale}}{$guess->{descid}};
    $guess2 = $h if (!defined $guess2);
    $tthen2 = $times{$guess2->{scale}};
    my $m = $guess->{orig_mvir};
    if ($guess2->{orig_mvir} != $m and $tthen != $tthen2) {
	$m += ($guess2->{orig_mvir} - $m)* (($tnow-$tthen)-$t) / ($tthen2-$tthen);
    }
    return ($m, $guess);
}

sub calc_accretion_rates {
    no warnings 'recursion';
    my $h = shift;
    return if ($h->{seen}>2);
    $h->{seen} = 3;
    if ($h->{prog}) {
	my ($h100, $hdyn, $h2dyn, $h3dyn, $h4dyn, $hmp) = calc_accretion_rates($h->{prog});
	$h4dyn = $h unless (defined($h4dyn));
	$h3dyn = $h unless (defined($h3dyn));
	$h2dyn = $h unless (defined($h2dyn));
	$hdyn = $h unless (defined($hdyn));
	$h100 = $h unless (defined($h100));
	$hmp = $h unless (defined($hmp));
	my $tdyn = $dyn_times{$h->{scale}};
	my $tzphalf = $zphalf_times{$h->{scale}};
	my ($mtdyn, $m2tdyn, $m3tdyn, $m4tdyn, $m100);
	($mtdyn, $hdyn) = find_mass_years_ago($h, $hdyn, $tdyn);
	($m2tdyn, $h2dyn) = find_mass_years_ago($h, $h2dyn, 2*$tdyn);
	($m3tdyn, $h3dyn) = find_mass_years_ago($h, $h3dyn, 3*$tdyn);
	($m4tdyn, $h4dyn) = find_mass_years_ago($h, $h4dyn, 4*$tdyn);
	($m100, $h100) = find_mass_years_ago($h, $h100, 100e6);
	(undef, $hmp) = find_mass_years_ago($h, $hmp, $tzphalf);
	$h->{acc_dyn} = ($h->{orig_mvir}-$mtdyn) / $tdyn;
	$h->{acc_2dyn} = ($h->{orig_mvir}-$m2tdyn) / (2*$tdyn);
	$h->{acc_3dyn} = ($h->{orig_mvir}-$m3tdyn) / (3*$tdyn);
	$h->{acc_4dyn} = ($h->{orig_mvir}-$m4tdyn) / (4*$tdyn);
	$h->{acc_100} = ($h->{orig_mvir}-$m100) / 100e6;
	$h->{acc_inst} = ($h->{orig_mvir} - $h->{prog}->{orig_mvir}) / ($times{$h->{scale}} - $times{$h->{prog}->{scale}});
	$h->{acc_mpeak} = ($h->{mpeak} - $hmp->{mpeak})/($times{$h->{scale}} - $times{$hmp->{scale}});
	$hdyn = $h unless (defined($hdyn));
	$h->{acc_dyn_dyn} = $hdyn->{acc_dyn};
	$h->{acc_2dyn_dyn} = $h2dyn->{acc_dyn};
	return ($h100, $hdyn, $h2dyn, $h3dyn, $h4dyn, $hmp);
    }
    else {
	my $tdyn = $dyn_times{$h->{scale}};
	$h->{acc_inst} = $h->{acc_100} = $h->{acc_dyn} = $h->{acc_2dyn} =  $h->{acc_3dyn} = $h->{acc_4dyn} = $h->{acc_dyn_dyn} = $h->{acc_2dyn_dyn} = $h->{acc_mpeak} = $h->{orig_mvir} / $tdyn;
	return ($h, $h, $h, $h, $h, $h);
    }
}


sub load_config {
    my $config = new ReadConfig($ARGV[0]);
    $config->set_defaults(
	SCALEFILE => "/Volumes/Peter 1/Bolshoi/DescScales.txt",
	TREE_OUTBASE => "/Volumes/Peter 2/Bolshoi/Trees",
	HLIST_OUTBASE => "/Volumes/Peter 2/Bolshoi/Hlists");
    our $TREE_OUTBASE = $config->{TREE_OUTBASE};
    our $HLIST_OUTBASE = $config->{HLIST_OUTBASE};
    our $OUTLIST = $config->{SCALEFILE};
    our $Om = $config->{Om} || 0.27;
    our $h0 = $config->{h0} || 0.70;
    our $MASS_RES_OK = $config->{MASS_RES_OK} || 1e11;
    Universe::Time::init($Om, $h0);
    
#    our $SINGLE_THREAD_OUTPUT = $config->{SINGLE_THREAD_OUTPUT};
}

package TreeHalo;

sub _open_scale {
    my $scale = shift;
    my $fn = $main::HLIST_OUTBASE.sprintf("/hlist_%.5f.list.%s", $scale,
	$main::filename);
    $tree_outputs{$scale} = undef;
    open $tree_outputs{$scale}, ">", $fn or
	die "Unable to open file $fn for writing!\n";
}

sub new {
    my $class = shift;
    my $obj = { map { $_ => 0 }
		qw(id descid scale desc_scale
num_prog pid upid desc_pid mtype mmp flags phantom
mvir rvir vmax vrms rs mtrunc orig_mvir rtrunc last_mm np spin) };
    $obj->{pos} = [0,0,0];
    $obj->{vel} = [0,0,0];
    $obj->{J} = [0,0,0];
    bless $obj, $class;
}

sub scanline {
    my $h = shift;
    my $line = shift;
    ($h->{scale}, $h->{id}, $h->{desc_scale}, $h->{descid}, $h->{num_prog},
     $h->{pid}, $h->{upid}, $h->{desc_pid}, $h->{phantom},
     $h->{mvir}, $h->{orig_mvir}, $h->{rvir}, $h->{rs}, $h->{vrms},
     $h->{mmp}, $h->{last_mm}, $h->{vmax},
     $h->{pos}[0], $h->{pos}[1], $h->{pos}[2],
     $h->{vel}[0], $h->{vel}[1], $h->{vel}[2],
     $h->{J}[0], $h->{J}[1], $h->{J}[2], $h->{spin}, $h->{rest}) = split(" ", $line, 28);
    chomp($h->{rest});
    $h->{mvir} = abs($h->{mvir});
    return "$h->{scale} $h->{id}";
}

sub print {
    my $h = shift;
    return unless (defined $h->{scale} and $h->{scale} > 0);
    _open_scale($h->{scale}) if (!exists($tree_outputs{$h->{scale}}));
    my $file = $tree_outputs{$h->{scale}};
    $file->printf("%.4f %8s %.4f %8s %6s %8s %8s %8s %2s %.5e %.5e %6f %6f %6f %2s %.4f %6f %.5f %.5f %.5f %.3f %.3f %.3f %.3e %.3e %.3e %.5f %s %.5e %.5e %6f %6f %.4f %.3e %.3e %.3e %.3e %.3e %.3e %.3f %.3f %.3e %.3f %.3f\n",
    $h->{scale}, $h->{id}, $h->{desc_scale}, $h->{descid}, $h->{num_prog},
    $h->{pid}, $h->{upid}, $h->{desc_pid}, $h->{phantom},
    $h->{mvir}, $h->{orig_mvir}, $h->{rvir}, $h->{rs}, $h->{vrms},
    $h->{mmp}, $h->{last_mm}, $h->{vmax},
    $h->{pos}[0], $h->{pos}[1], $h->{pos}[2],
    $h->{vel}[0], $h->{vel}[1], $h->{vel}[2],
    $h->{J}[0], $h->{J}[1], $h->{J}[2], $h->{spin}, $h->{rest},
    $h->{macc}, $h->{mpeak}, $h->{vacc}, $h->{vpeak}, $h->{halfmass},
    $h->{acc_inst}, $h->{acc_100}, $h->{acc_dyn}, $h->{acc_2dyn}, $h->{acc_mpeak}, $h->{mpeak_scale}, $h->{acc_scale}, $h->{first_acc}{scale}, $h->{first_acc}{orig_mvir}, $h->{first_acc}{vmax}, $h->{vmpeak});
}

