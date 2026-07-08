
STDOUT.sync = true; STDERR.sync = true

require 'fileutils'
require 'shellwords'

codon_table, fque, faa, fout, flog = ARGV

raise("usage: ruby predict_nucl_query.rb <codon_table> <query nucleotide fasta> <protein fasta> <gff output> <log>") unless flog

FileUtils.mkdir_p(File.dirname(faa))

tmp_input = nil
input = fque

if fque =~ /\.gz$/ or fque =~ /\.gzip$/
  require 'zlib'
  tmp_input = File.join(File.dirname(faa), "#{File.basename(faa, ".faa")}.prodigal_input.fna")
  Zlib::GzipReader.open(fque) do |fr|
    open(tmp_input, "w") do |fw|
      while chunk = fr.read(1024 * 1024)
        fw.write(chunk)
      end
    end
  end
  input = tmp_input
end

cmd = ["prodigal", "-p", "meta", "-g", codon_table, "-i", input, "-a", faa, "-o", fout]

open(flog, "w") do |fw|
  fw.puts "$ #{cmd.shelljoin}"
  fw.flush
  ok = system(*cmd, out: fw, err: fw)
  raise("prodigal failed with status #{$?.exitstatus}. See #{flog}") unless ok
end

FileUtils.rm_f(tmp_input) if tmp_input
