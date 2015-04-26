package MultiThread;
use IO::Select;

our @pids;
our $threadlimit = 16;

#$SIG{CHLD} = \&WaitForChild;

sub ThreadLimit {
    my $new = shift;
    $threadlimit = $new if ($new);
    return $threadlimit;
}

sub EvalChildren {
    my @execs = @_;
    #To be completed...
}

sub ForkChild {
    WaitForChild() until (@pids < $threadlimit);
    my $pid = fork();
    if (!$pid) {
	@pids = ();
	return 1;
    }
    push @pids, $pid;
    return;
}

sub ExecChild {
    return unless ForkChild();
    exec(@_);
    exit;
}

sub ClearPids {
    @pids = ();
}

sub WaitForChild {
    return unless @pids;
    my $pid = wait();
    @pids = grep { $_ != $pid } @pids;
}

sub WaitForChildren {
    while (@pids) {
	my $pid = shift @pids;
	waitpid($pid, 0);
    }
}


1;
