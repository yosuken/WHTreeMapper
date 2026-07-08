
STDOUT.sync = true; STDERR.sync = true

# {{{ procedures
WriteBatch  = lambda do |t, jobdir, outs|
  jdir = "#{jobdir}/#{t.name.split(":")[-1]}"; mkdir_p jdir unless File.directory?(jdir)
  jnum = outs.size

  if jnum > 0
    outs.each_slice(jnum).with_index(1){ |ls, idx| ## always 1 file
      open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.sh", "w"){ |fjob|
        fjob.puts ls
      }
    }
  end
end

RunBatch  = lambda do |t, jobdir, ncpu, logdir|
  jdir = "#{jobdir}/#{t.name.split(":")[-1]}"
  ldir = "#{logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)
  tdir = "#{jdir}/tmp"; mkdir_p tdir unless File.directory?(tdir)

  Dir["#{jdir}/*.sh"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin| ## always 1 or 0 file
    sh "TMPDIR=#{tdir} parallel --jobs #{ncpu} --joblog #{ldir}/parallel.log <#{fin}"
  }
  open("#{ldir}/exit", "w"){ |fw| fw.puts "exit with status 0 at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
end

PrintStatus = lambda do |current, total, status, t|
  puts ""
  puts "\e[1;32m===== #{Time.now}\e[0m"
  puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
  puts ""
end

CheckVersion = lambda do |commands|
  commands.each{ |command|
    str = case command
    when "ruby"
      %{which ruby && ruby -v 2>&1}
    when "hmmsearch"
      %{which hmmsearch && hmmsearch -h 2>&1 |head -n 2}
    when "parallel"
      %{which parallel && LANG=C parallel --version 2>&1 |head -n 1}
    when "diamond"
      %{which diamond && diamond version 2>&1}
    when "prodigal"
      %{which prodigal && prodigal -v 2>&1 |head -n 2}
    end
    puts ""
    puts "\e[1;32m===== check version: #{command}\e[0m"
    puts ""
    puts "$ #{str}"
    puts `#{str}`
  }
end
# }}} procedures


# {{{ task controller
task :default do
  ### constants
  Errmsg    = "\e[1;31mError:\e[0m"
  Warmsg    = "\e[1;35mWarning:\e[0m"
  Sccmsg    = "\e[1;32mSuccess:\e[0m"

  ### constants from arguments
  Odir         = ENV["outdir"]
  OdirExist    = ENV["outdir_exist"]
  Fques        = ENV["input"]
  InputMode    = ENV["input_mode"]
  raise("#{Errmsg} input_mode should be prot or nucl.") unless %w|prot nucl|.include?(InputMode)
  CodonTable   = (ENV["codon_table"] || "11").to_i
  Rpkgs        = ENV["refpkg"]
  Ncpu         = ENV["ncpus"].to_i
  Evalue       = ENV["evalue"].to_f
  EvalueDom    = ENV["evaluedom"].to_f

  MinSeqLen    = ENV["minseqlen"].to_i
  MinHmmLen    = ENV["minhmmlen"].to_i
  MinHmmCov    = ENV["minhmmcov"].to_f
  MinAliLen    = ENV["minalilen"].to_i
  MinAliCov    = ENV["minalicov"].to_f
  MinHmmLenDom = ENV["minhmmlendom"].to_i
  MinHmmCovDom = ENV["minhmmcovdom"].to_f
  MinAliLenDom = ENV["minalilendom"].to_i
  MinAliCovDom = ENV["minalicovdom"].to_f

  DiamondDB    = ENV["diamond_db"]
  DmndEvalue   = ENV["dmnd_evalue"].to_f
  DmndId       = ENV["dmnd_id"].to_f
  DmndScov     = ENV["dmnd_subject_cover"].to_f
  DmndMaxTgt   = ENV["dmnd_max_target_seqs"].to_i
  Z            = 1_000_000 ## hmmsearch data size for evalue calculation

  ### define tasks
  tasks = []
  tasks << "01-1a-A.validate_query"
  tasks << "01-1a-B.parse_query_info"
  tasks << "01-1b-A.validate_refpkg"
  tasks << "01-1b-B.parse_refpkg_info"
  tasks << "01-2a.hmmsearch"
  tasks << "01-2b.parse_hmmsearch"
  tasks << "01-2c.split_fasta"
  tasks << "01-2d.copy_detected"
  tasks << "02-1.merge_regions"
  tasks << "02-2.diamond_blastp"
  tasks << "02-3.create_detected_tsv"

  ### store Tasks
  Tasks = tasks

  ### check version
  commands = %w|hmmsearch ruby parallel diamond|
  commands << "prodigal" if InputMode == "nucl"
  CheckVersion.call(commands)

  ## dir path
  Jobdir    = "#{Odir}/batch"
  Pkgdir    = "#{Odir}/refpkg"
  Predir    = "#{Odir}/prefilter"
  PreQuedir = "#{Odir}/prefilter/query"
  PreFildir = "#{Odir}/prefilter/hmmsearch"
  Resdir    = "#{Odir}/hmm_hits"
  Dmnddir   = "#{Odir}/diamond"
  Logdir    = "#{Odir}/log/tasks"
  SeqInfo   = ENV["seq_info"]

  ## variables
  $fques    = []
  $qnames   = {}
  $fpkgs    = []
  $rnames   = {}

  ## Odir exist? (warning is already printed by WHTreeMapper entry point)

  ### run
  NumStep  = Tasks.size
  Tasks.each.with_index(1){ |task, idx|
    Rake::Task[task].invoke(idx)
  }
end
# }}} default (run all tasks)


# {{{ tasks
# {{{ desc "01-1a-A.validate_query"
desc "01-1a-A.validate_query"
task "01-1a-A.validate_query", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs          = []
  script        = "#{__dir__}/script/#{t.name}.rb"
  script_nucl   = "#{__dir__}/script/predict_nucl_query.rb"
  query_name_of = lambda do |fque|
    name = File.basename(fque).gsub(/\.gz$/, "").gsub(/\.gzip$/, "").split(".")[0..-2]*"."
    name.empty? ? File.basename(fque).gsub(/\.gz$/, "").gsub(/\.gzip$/, "") : name
  end

  fques = Fques.split(/[,\s]+/).inject([]){ |a, path| a += Dir[path.gsub("~", ENV["HOME"])].sort }
  fques = fques.map{ |fque|
    if File.zero?(fque)
      $stderr.puts "#{Warmsg} #{fque} is empty. Skip it."
      next nil
    else
      File.open(fque){ |fr|
        if fque =~ /\.gz$/ or fque =~ /\.gzip/
          is_fque_gz = true
          require 'zlib'

          Zlib::GzipReader.open(fque) { |fr|
            l = fr.gets
            raise("#{Errmsg} <query> file: #{fque} has wrong format.\nThe first line should begin '>'") if l[0] != ">"
          }
        else
          l = fr.gets
          raise("#{Errmsg} <query> file: #{fque} has wrong format.\nThe first line should begin '>'") if l[0] != ">"
        end
      }
    end
    fque
  }.compact

  $stderr.puts ["", "", "\e[1;32m===== check query file (N=#{fques.size}) \e[0m"]
  raise("#{Errmsg} no query file detected.") if fques.size == 0

  mkdir_p PreQuedir

  if InputMode == "nucl"
    pdir = "#{PreQuedir}/prodigal"; mkdir_p pdir
    prodigal_outs = []
    fques = fques.map{ |fque|
      name = query_name_of.call(fque)
      faa  = "#{pdir}/#{name}.faa"
      fout = "#{pdir}/#{name}.gff"
      flog = "#{pdir}/#{name}.prodigal.log"
      prodigal_outs << "ruby #{script_nucl} #{CodonTable} #{fque} #{faa} #{fout} #{flog}" unless File.exist?(faa)
      faa
    }

    if prodigal_outs.size > 0
      WriteBatch.call(t, Jobdir, prodigal_outs)
      RunBatch.call(t, Jobdir, Ncpu, Logdir)
    end

    fques = fques.map{ |fque|
      if File.zero?(fque)
        $stderr.puts "#{Warmsg} #{fque} is empty. Skip it."
        next nil
      end
      fque
    }.compact
    raise("#{Errmsg} no protein sequence was predicted from nucleotide query file.") if fques.size == 0
  end

  fques.each{ |fque|
    name = query_name_of.call(fque)
    raise("#{Errmsg} file name #{name} is given twice.") if $qnames[name]
    $qnames[name] = 1
  }

  flsts  = []
  fques.each{ |fque|
    name = query_name_of.call(fque)
    fa   = "#{PreQuedir}/#{name}.fa"
    fjsn = "#{PreQuedir}/#{name}.json"
    flst = "#{PreQuedir}/#{name}.list"

    if !File.exist?(fjsn)
      outs << "ruby #{script} #{MinSeqLen} #{name} #{fque} #{fa} #{fjsn} #{flst}"
    end

    flsts  << flst
  }

  next if outs.size == 0

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)

  sh "cat #{flsts*' '} |sort |uniq -d >#{Predir}/duplicated.list"
  n_dup = IO.readlines("#{Predir}/duplicated.list").size

  puts "#{Warmsg} #{n_dup} non-unique sequence IDs are found between query files. See #{Predir}/duplicated.list" if n_dup > 0
end
# }}}

# {{{ desc "01-1a-B.parse_query_info"
desc "01-1a-B.parse_query_info"
task "01-1a-B.parse_query_info", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  require 'json'

  $fques = []
  $qnames.each_key{ |name|
    fjsn = "#{PreQuedir}/#{name}.json"
    sjsn = IO.readlines(fjsn)[0]
    puts sjsn

    $fques << JSON.parse(sjsn, symbolize_names: true)
  }
end
# }}}

# {{{ desc "01-1b-A.validate_refpkg"
desc "01-1b-A.validate_refpkg"
task "01-1b-A.validate_refpkg", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []
  script = "#{__dir__}/script/01-1b-A.validate_refpkg_detect.rb"

  rpkgs = Rpkgs.split(/[,\s]+/).sort_by{ |path| File.basename(path) }.inject([]){ |a, path| a += Dir[path.gsub("~", ENV["HOME"])].sort }
  $stderr.puts ["", "", "\e[1;32m===== check refpkg (N=#{rpkgs.size}) \e[0m"]
  raise("#{Errmsg} no refpkg directory detected.") if rpkgs.size == 0

  $rnames = {}
  rpkgs.each{ |rpkg|
    name = File.basename(rpkg)
    raise("#{Errmsg} refpkg name #{name} is not unique.") if $rnames[name]
    $rnames[name] = 1

    fa    = Dir["#{rpkg}/*.fa"] + Dir["#{rpkg}/*.mfa"] + Dir["#{rpkg}/*.fasta"] + Dir["#{rpkg}/*.faa"]
    fhmm  = Dir["#{rpkg}/*.hmm"]
    ftre  = Dir["#{rpkg}/*.tree"] + Dir["#{rpkg}/*.nwk"] + Dir["#{rpkg}/*.newick"]

    raise("#{Errmsg} #{rpkg} is not a directory.") unless File.directory?(rpkg)
    raise("#{Errmsg} #{rpkg} does not contain fasta file. #{rpkg}/*{.fa|.mfa|.faa|.fasta} should exist.") if fa.size == 0
    raise("#{Errmsg} #{rpkg} contains multiple fasta files. #{rpkg}/*{.fa|.mfa|.faa|.fasta} should be only one.") if fa.size > 1
    raise("#{Errmsg} #{rpkg} contains multiple tree files. #{rpkg}/*{.tree|.nwk|.newick} should be only one.") if ftre.size > 1
    raise("#{Errmsg} #{rpkg} does not contain .hmm file. #{rpkg}/*.hmm should exist.") if fhmm.size == 0
    raise("#{Errmsg} #{rpkg} contains multiple .hmm files. #{rpkg}/*.hmm should be only one.") if fhmm.size > 1

    falnO  = fa[0]
    ftreO  = ftre[0]
    fhmmO  = fhmm[0]

    odir = "#{Pkgdir}/#{name}"; mkdir_p odir unless File.directory?(odir)
    faln = "#{Pkgdir}/#{name}/backbone.mfa"
    flog = "#{Pkgdir}/#{name}/backbone.log"
    if !File.exist?(faln)
      outs << "ruby #{script} #{rpkg} #{odir} #{falnO} #{ftreO} #{fhmmO} >#{flog} 2>&1"
    end
  }

  next if outs.size == 0

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-1b-B.parse_refpkg_info"
desc "01-1b-B.parse_refpkg_info"
task "01-1b-B.parse_refpkg_info", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  require 'json'

  $rnames.each_key{ |name|
    fjsn = "#{Pkgdir}/#{name}/backbone.json"
    sjsn = IO.readlines(fjsn)[0]
    puts sjsn

    $fpkgs << JSON.parse(sjsn, symbolize_names: true)
  }
end
# }}}

# {{{ desc "01-2a.hmmsearch"
desc "01-2a.hmmsearch"
task "01-2a.hmmsearch", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  outs   = []

  npara = [Ncpu, $fpkgs.size].min
  ncpu  = [1, (Ncpu.to_f / npara).round].max

  puts "###"
  puts "### number of parallel for hmmsearch: #{npara}"
  puts "### number of cpu per a hmmsearch: #{ncpu} CPUs (total: #{Ncpu})"
  puts "###"

  $fques.each{ |que|
    next if que[:numseq] == 0
    odir = "#{PreFildir}/#{que[:name]}/out"; mkdir_p odir
    $fpkgs.each{ |pkg|
      db   = pkg[:fhmm]
      fa   = que[:fasta]
      fout = "#{odir}/#{pkg[:name]}.out"
      flog = "#{odir}/#{pkg[:name]}.log"
      outs << "hmmsearch --cpu #{ncpu} -Z #{Z} --notextw -o #{fout} #{db} #{fa} >#{flog} 2>&1"
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, npara, Logdir)
end
# }}}

# {{{ desc "01-2b.parse_hmmsearch"
desc "01-2b.parse_hmmsearch"
task "01-2b.parse_hmmsearch", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  outs   = []

  script = "#{__dir__}/script/parse_hmmsearch.rb"

  $fques.each{ |que|
    fa    = que[:fasta]
    idir  = "#{PreFildir}/#{que[:name]}/out"
    next if Dir["#{idir}/*.out"].size == 0

    odir = "#{PreFildir}/#{que[:name]}/parsed"; mkdir_p odir unless File.directory?(odir)
    fhmm = "#{odir}/hmmsearch_concat.out"
    open(fhmm, "w"){ |fw|
      $fpkgs.each{ |pkg|
        f = "#{idir}/#{pkg[:name]}.out"
        if File.exist?(f) and !File.zero?(f)
          fw.puts IO.read(f)
        else
          raise("#{Errmsg} #{f} does not exist or zero-sized.")
        end
      }
    }

    flog    = "#{odir}/parse.log"
    option  = "-ge #{Evalue} -e #{EvalueDom} --create-evalue-table --min-hmm-len-dom #{MinHmmLenDom} --min-hmm-cov-dom #{MinHmmCovDom} --min-ali-len-dom #{MinAliLenDom} --min-ali-cov-dom #{MinAliCovDom} "
    option += "--min-hmm-len #{MinHmmLen} --min-hmm-cov #{MinHmmCov} --min-ali-len #{MinAliLen} --min-ali-cov #{MinAliCov}"
    outs << "ruby #{script} #{option} -i #{fhmm} -f #{fa} -o #{odir} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-2c.split_fasta"
desc "01-2c.split_fasta"
task "01-2c.split_fasta", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")

  # {{{ def parse_fasta(fin, fout, ids)
  def parse_fasta(fin, fout, ids)
    open(fout, "w"){ |fw|
      flag, ent = 0, []
      open(fin){ |fr|
        while l = fr.gets
          if l[0] == ">"
            fw.puts ent if ent != []
            flag, ent = 0, []

            id = l.strip[1..-1].split(/\s+/)[0]
            if ids[id]
              flag = 1
              ent << l.strip
            end
          else
            ent << l.strip if flag == 1
          end
        end
        fw.puts ent if ent != []
      }
    }
  end
  # }}}

  $fques.each{ |que|
    ftsv = "#{PreFildir}/#{que[:name]}/parsed/best-hit.tsv"
    hmm2regs = {}
    hmm2gids = {}
    open(ftsv){ |fr|
      _ = fr.gets # skip header
      while l = fr.gets
        a = l.strip.split("\t", -1)
        gid, hmm, reg = a.values_at(0, 3, 21)
        hmm2gids[hmm] ||= {}
        hmm2gids[hmm][gid] = 1
        hmm2regs[hmm] ||= {}
        hmm2regs[hmm][reg] = 1
      end
    }

    $fpkgs.each{ |pkg|
      regs = hmm2regs[pkg[:hmmname]]
      gids = hmm2gids[pkg[:hmmname]]
      next if !regs or regs.keys.size == 0

      %w|whole region|.zip(%w|best-hit.whole.fa best-hit.fa|, [gids, regs]){ |type, fname, ids|
        fin  = "#{PreFildir}/#{que[:name]}/parsed/#{fname}"
        odir = "#{PreFildir}/#{que[:name]}/seq/#{type}/#{pkg[:name]}"; mkdir_p odir unless File.directory?(odir)
        fout = "#{odir}/#{que[:name]}.fa"

        parse_fasta(fin, fout, ids)
      }
    }
  }
end
# }}}

# {{{ desc "01-2d.copy_detected"
desc "01-2d.copy_detected"
task "01-2d.copy_detected", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  outs = []

  $fpkgs.each{ |pkg|
    %w|whole region|.each{ |type|
      fins = []
      $fques.each{ |que|
        fin  = "#{PreFildir}/#{que[:name]}/seq/#{type}/#{pkg[:name]}/#{que[:name]}.fa"
        next unless File.exist?(fin)
        fins << fin
      }

      next if fins.size == 0
      odir = "#{Resdir}/#{pkg[:name]}/seq"; mkdir_p odir unless File.directory?(odir)
      fout = "#{odir}/#{type}.fa"
      outs << "cat #{fins*' '} >#{fout}"
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "02-1.merge_regions"
desc "02-1.merge_regions"
task "02-1.merge_regions", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")

  odir = "#{Odir}/diamond/query"
  mkdir_p odir unless File.directory?(odir)
  fout = "#{odir}/region.faa"

  nseq = 0
  open(fout, "w") { |fw|
    $fpkgs.each { |pkg|
      fin = "#{Resdir}/#{pkg[:name]}/seq/region.fa"
      next unless File.exist?(fin) && !File.zero?(fin)

      clade = pkg[:name]
      IO.read(fin).split(/^>/)[1..-1].each { |ent|
        lab, *seq = ent.split("\n")
        gid = lab.split(" ")[0]
        fw.puts ">#{gid} #{clade}"
        fw.puts seq.join("").strip
        nseq += 1
      }
    }
  }

  if nseq == 0
    $stderr.puts "#{Warmsg} No WH regions detected in any clade. Diamond blastp will be skipped."
  else
    puts "### Merged #{nseq} region sequences into #{fout}"
  end

  ldir = "#{Logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)
  open("#{ldir}/exit", "w"){ |fw| fw.puts "exit with status 0 at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
end
# }}}

# {{{ desc "02-2.diamond_blastp"
desc "02-2.diamond_blastp"
task "02-2.diamond_blastp", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")

  fquery = "#{Odir}/diamond/query/region.faa"

  unless File.exist?(fquery) && !File.zero?(fquery)
    puts "#{Warmsg} No query file for diamond blastp. Skipped."
    ldir = "#{Logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)
    open("#{ldir}/exit", "w"){ |fw| fw.puts "exit with status 0 (skipped) at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
    next
  end

  mkdir_p Dmnddir unless File.directory?(Dmnddir)
  fout = "#{Dmnddir}/dmnd.blastp.out"
  flog = "#{Dmnddir}/dmnd.blastp.log"

  outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp"

  outs = []
  outs << "diamond blastp --ultra-sensitive --id #{DmndId} --subject-cover #{DmndScov} --threads #{Ncpu} --max-target-seqs #{DmndMaxTgt} --dbsize 1e9 --evalue #{DmndEvalue} --outfmt #{outfmt} --query #{fquery} --db #{DiamondDB} --out #{fout} >#{flog} 2>&1"

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, 1, Logdir)
end
# }}}

# {{{ desc "02-3.create_detected_tsv"
desc "02-3.create_detected_tsv"
task "02-3.create_detected_tsv", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")

  fdmnd = "#{Dmnddir}/dmnd.blastp.out"

  unless File.exist?(fdmnd) && !File.zero?(fdmnd)
    puts "#{Warmsg} No diamond blastp result. Skipped."
    ldir = "#{Logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)
    open("#{ldir}/exit", "w"){ |fw| fw.puts "exit with status 0 (skipped) at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
    next
  end

  ### load seq_info mapping: region => {clade, subclade, index_in_tree, leaf_name}
  ### build both exact match and GID-prefix match (strip _fmX_toY)
  region2info = {}
  gid2info    = {}
  open(SeqInfo){ |fr|
    fr.gets # skip header
    while l = fr.gets
      a = l.chomp.split("\t", -1)
      clade, subclade, idx, leaf_name, region = a
      info = { clade: clade, subclade: subclade, index_in_tree: idx, leaf_name: leaf_name }
      region2info[region] = info
      gid = region.sub(/_fm\d+_to\d+$/, "")
      gid2info[gid] ||= info
    end
  }

  fout = "#{Odir}/detected.tsv"

  nout = 0
  open(fout, "w"){ |fw|
    fw.puts %w|query detected_wh_region clade subclade index_in_tree leaf_name qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp|.join("\t")

    open(fdmnd){ |fr|
      while l = fr.gets
        a = l.chomp.split("\t", -1)
        qseqid = a[0]
        sseqid = a[1]

        ### qseqid has the detected region suffix; query is the original protein ID.
        query = qseqid.sub(/_fm\d+_to\d+$/, "")

        ### extract detected WH region from last _fmX_toY in qseqid
        wh_region = (qseqid =~ /_fm(\d+)_to(\d+)$/) ? "#{$1}-#{$2}" : ""

        info = region2info[sseqid]
        if !info
          gid = sseqid.sub(/_fm\d+_to\d+$/, "")
          info = gid2info[gid]
        end
        if info
          fw.puts [query, wh_region, info[:clade], info[:subclade], info[:index_in_tree], info[:leaf_name], a].flatten.join("\t")
        else
          fw.puts [query, wh_region, "", "", "", "", a].flatten.join("\t")
          $stderr.puts "#{Warmsg} sseqid #{sseqid} not found in seq_info."
        end
        nout += 1
      end
    }
  }

  puts "### Created #{fout} with #{nout} entries"

  ldir = "#{Logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)
  open("#{ldir}/exit", "w"){ |fw| fw.puts "exit with status 0 at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
end
# }}}

# }}} tasks
