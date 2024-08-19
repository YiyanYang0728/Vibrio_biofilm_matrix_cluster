#!/usr/bin/perl
$host = shift;
$instance = shift;
$arg = shift;

#### random sleep, rand() can be a fraction of second
select(undef,undef,undef,rand());

if ($arg) {
  @ids = split(/,/, $arg);
}
else {
  while(1) {
    if (opendir(DDIR, "all_pro.faa_40-seq")) { 
      @ids = grep {/^\d+$/} readdir(DDIR);
      last;
    }
    else {
      sleep(1);
    }
  }
}

foreach $id (@ids) {

  next unless (-e "all_pro.faa_40-seq/$id");
  next if (-e "all_pro.faa_40-seq/$id.lock");
  $cmd = `touch all_pro.faa_40-seq/$id.lock`;

  if (50) {
    $cmd = `blastp -outfmt 6 -db ./all_pro.faa_40.1354641 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query all_pro.faa_40-seq/$id -out all_pro.faa_40-bl/$id`;
    $cmd =                         `/usr/local/apps/cd-hit/cdhit-4.8.1/psi-cd-hit/psi-cd-hit.pl -J parse_blout_multi all_pro.faa_40-bl/$id -c .30000000000000000000 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 0 >> all_pro.faa_40-blm/$host.$instance`;
  }
  elsif (1) {
    $cmd = `blastp -outfmt 6 -db ./all_pro.faa_40.1354641 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query all_pro.faa_40-seq/$id | /usr/local/apps/cd-hit/cdhit-4.8.1/psi-cd-hit/psi-cd-hit.pl -J parse_blout all_pro.faa_40-bl/$id -c .30000000000000000000 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 1`;
  }
  else {
    $cmd = `blastp -outfmt 6 -db ./all_pro.faa_40.1354641 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query all_pro.faa_40-seq/$id -out all_pro.faa_40-bl/$id`;
    $cmd =                         `/usr/local/apps/cd-hit/cdhit-4.8.1/psi-cd-hit/psi-cd-hit.pl -J parse_blout all_pro.faa_40-bl/$id -c .30000000000000000000 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 0`;
  }
  $cmd = `rm -f  all_pro.faa_40-seq/$id`;
  $cmd = `rm -f  all_pro.faa_40-seq/$id.lock`;
}

($tu, $ts, $cu, $cs) = times();
$tt = $tu + $ts + $cu + $cs;
$cmd = `echo $tt >> all_pro.faa_40-seq/host.$host.$instance.cpu`;

