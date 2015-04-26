#!/usr/bin/perl -w
package CachedReadWrite;
use Symbol;
use Fcntl qw(:flock);
use IO::Socket::INET;
use IO::Handle;

sub new {
    my ($class, %opts) = @_;
    my $obj = gensym;
    %{*$obj} = (cachesize => 10000, rpos => 0,
		wpos => 0, wcache => [],
		%opts);
    @{*$obj} = (); #rcache
    bless $obj, $class;
    tie *{$obj}, "CachedReadWrite::Handle", $obj;
    return $obj;
}

sub set_append {    
    my $obj = shift;
    my $hash = \%{*$obj};
    $hash->{append} = shift;
}

sub _flush_wcache {
    my $obj = shift;
    my $hash = \%{*$obj};
    if (@{$hash->{wcache}}) {
	local *FILE;
	if ($hash->{append}) {
	    open FILE, ">>", $hash->{file} or
		die ("Couldn't open file $hash->{file} for appending!\n");
	} else {
	    unless (open FILE, "+<", $hash->{file}) {
		open FILE, ">", $hash->{file} or
		    die ("Couldn't open file $hash->{file} for writing!\n");
	    }
	    seek FILE, $hash->{wpos}, 0;
	}
	flock FILE, LOCK_EX;
	if ($hash->{append}) { seek FILE, 0, 2; }
	if ((-z FILE) and defined($hash->{header})) {
	    print FILE $hash->{header};
	}
	print FILE @{$hash->{wcache}};
	flock FILE, LOCK_UN;
	close FILE;
	$hash->{wpos} += length($_) for (@{$hash->{wcache}});
	@{$hash->{wcache}} = ();
    }
}

sub _fill_rcache {
    my $obj = shift;
    my $hash = \%{*$obj};
    my $cache = \@{*$obj};
    local (*FILE, $_);
    open FILE, "<", $hash->{file} or
	die ("Couldn't open file $hash->{file} for reading!\n");
    seek FILE, $hash->{rpos}, 0;
    while (defined($_=<FILE>)) {
	push @$cache, $_;
	$hash->{rpos} += length($_);
	last if (@$cache >= $hash->{cachesize});
    }
    close FILE;
}

sub rpos {
    my $obj = shift;
    my $hash = \%{*$obj};
    $hash->{rpos};
}

sub close {
    my $self = shift;
    $self->_flush_wcache();
    local @_;
    untie *$self;
}

sub print {
    my $self = shift;
    my $hash = \%{*$self};
    push @{$hash->{wcache}}, @_;
    $self->_flush_wcache() if (@{$hash->{wcache}} > $hash->{cachesize});
    return 1;
}

sub printf {
    my ($self,$format) = (shift,shift);
    local $\;
    return $self->print(sprintf($format, @_));
}

sub getc {
    my $self = shift;
    my $hash =  \%{*$self};
    my $cache = \@{*$self};
    $self->_fill_rcache() unless (@$cache);
    return unless (@$cache);
    return unless (length($cache->[0]));
    my $c = substr($cache->[0], 0, 1, "");
    shift @$cache unless (length($cache->[0]));
    return $c;
}

sub read {
    die "Unimplemented. :(\n";
}

sub readline {
    my $self = shift;
    my $cache = \@{*$self};
    $self->_fill_rcache() unless (@$cache);
    return unless (@$cache);
    return shift(@$cache);
}

sub unreadline {
    my $self = shift;
    #my $cache = \@{*$self};
    return unshift(@{*$self}, $_[0]);
}

sub DESTROY {
    my $self = shift;
    $self->_flush_wcache();
}

package CachedReadWrite::Handle;
use Scalar::Util qw(weaken);

sub TIEHANDLE {
    my ($class, $handle) = @_;
    weaken($handle);
    bless \$handle, $class;
}


sub READ     { ${shift()}->read     (@_) }
sub READLINE { ${shift()}->readline (@_) }
sub GETC     { ${shift()}->getc     (@_) }

sub PRINT    { ${shift()}->print    (@_) }
sub PRINTF   { ${shift()}->printf   (@_) }
sub WRITE    { ${shift()}->write    (@_) }

sub FILENO   { ${shift()}->fileno   (@_) }

sub CLOSE {                          #<---- Do not change this function!
    my $obj = ${$_[0]};
    local @_;
    $obj->close();
}


1;
