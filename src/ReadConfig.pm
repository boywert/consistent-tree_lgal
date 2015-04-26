#!/usr/bin/perl -w
package ReadConfig;

sub new {
    my ($class, $file) = @_;
    return bless {}, $class unless defined $file;
    open FILE, "<", $file or 
	die "Couldn't open file $file for reading!\n";
    local $_;
    my %config;
    while (<FILE>) {
	s/\#.*//;
	my ($key, $value) = /([^=]+)=(.*)/;
	next unless defined $value;
	print "$key = $value\n";
	($key, $value) = map { trim($_) } ($key, $value);
	next unless length($key) and length($value);
	$config{$key} = $value;
    }
    close FILE;
    bless \%config, $class;
}

sub trim {
    my $val = shift;
    $val =~ s/^\s*['"]?//;
    $val =~ s/['"]?\s*$//;
    return $val;
}

sub set_defaults {
    my $self = shift;
    my %opts = @_;
    %$self = (%opts, %$self);
}

sub print_config {
    my ($self, $fh) = @_;
    $fh = \*STDOUT unless $fh;
    for (sort keys %$self) {
	print $fh "\"$_\" = \"$self->{$_}\"\n";
    }
}

1;
