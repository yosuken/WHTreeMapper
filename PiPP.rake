

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

  Dir["#{jdir}/*.sh"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin| ## always 1 or 0 file
    if ncpu > 1
      sh "parallel --jobs #{ncpu} --joblog #{ldir}/parallel.log <#{fin}"
    else
      # sh "bash -c #{fin}" ## --> permission denied
      sh "bash #{fin}"
    end
  }
  open("#{ldir}/exit", "w"){ |fw| fw.puts "exit at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
end

PrintStatus = lambda do |current, total, status, t|
  puts ""
  puts "\e[1;32m===== #{Time.now}\e[0m"
  puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
  puts ""
  $stdout.flush
end

CheckVersion = lambda do |commands|
  commands.each{ |command|
    str = case command
    when "ruby"
      %|which ruby && ruby --version 2>&1|
    when "hmmsearch"
      %{which hmmsearch && hmmsearch -h 2>&1 |head -n 2}
    when "mafft"
      %{which mafft && mafft -v 2>&1 |head -n 4 |tail -n 1}
    when "gappa"
      %{which gappa && gappa --version}
    when "pplacer"
      %{which pplacer && pplacer --version}
    when "parallel"
      %{which parallel && LANG=C parallel --version 2>&1 |head -n 1}
    end
    puts ""
    puts "\e[1;32m===== check version: #{command}\e[0m"
    puts ""
    puts "$ #{str}"
    ### run
    puts `#{str}`
    ### flush
    $stdout.flush
  }
end
# }}} procedures


# {{{ task controller
task :default do
  ### define tasks
  tasks = []
  tasks << "01-1a-A.validate_query"
  tasks << "01-1a-B.parse_query_info"
  tasks << "01-1b.validate_refpkg"
  tasks << "01-2a.hmmsearch"
  tasks << "01-2b.parse_hmmsearch"
  tasks << "01-2c.merge_fasta"
  tasks << "01-3a.chunkify"
  tasks << "01-3b.mafft_add"
  tasks << "01-3c.pplacer"
  tasks << "01-3d.unchunkify"
  tasks << "01-3e.unchunkify_alignment"
  tasks << "01-4a.info"
  tasks << "01-4b.lwr"
  tasks << "01-4d.assign"  ## if ftax is given
  tasks << "01-4f.graft"
  tasks << "01-4h.aligned_position" ## if fpos is given
  tasks << "01-6a.aa_feature"

  # tasks << "01-4e.extract" ## if ftax is given ### memory requirement is too high for rep_2022-06-08
  # tasks << "01-4g.heat-tree"
  # tasks << "01-5a.krd"

  ### currently not used
  ### tasks << "01-4c.edpl"    ## high memory requirement. --> do not run
  ### tasks << "01-5b.edgepca" ## high memory requirement if input (num of samples?) is large.
  ### tasks << "01-5c.squash"
  ### tasks << "01-5d.dispersion"

  if ENV["only_detect"] == "true"
    Tasks = tasks[0, 5]
  else
    Tasks = tasks
  end

  ### constants from arguments
  Odir         = ENV["outdir"]            ## output directory
  OdirExist    = ENV["outdir_exist"]      ## output directory exist? ("true" or "")
  Fques        = ENV["query"]             ## query
  Rpkgs        = ENV["refpkg"]            ## refpkg
  Ncpu         = ENV["ncpus"].to_i        ## num of CPUs
  Evalue       = ENV["evalue"].to_f       ## hmmsearch evalue threshold (default: 1e-5)
  MinSeqLen    = ENV["minseqlen"].to_i    ## minimum query length
  MinHmmLen    = ENV["minhmmlen"].to_i    ## hmmsearch minimum hit hmm length
  MinHmmLenFrc = ENV["minhmmlenfrc"].to_f ## minimum fraction of hmmsearch hit hmm length
  C_size       = ENV["chunk_size"].to_i   ## chunk size
  M_mafft      = ENV["mafft_method"]      ## FFT-NS-2 or FFT-NS-i or E-INS-i
  Ex_lvs       = ENV["extract_levels"]    ## clade/taxonomy level for 'gappa prepare extract' (default: 0)
  Trim_opt     = ENV["trim_option"]       ## 'merge' (merge regions of hmmserach hits)  or 'largest' (take largest hit) (default: merge)

  ### check version
  commands  = %w|hmmsearch mafft gappa pplacer ruby|
  commands += %w|parallel| if Ncpu > 1
  CheckVersion.call(commands)

  ### constants
  Z         = 1_000_000 # hmmsearch data size for evalue calculation
  Errmsg    = "\e[1;31mError:\e[0m"
  Warmsg    = "\e[1;35mWarning:\e[0m"

  ## dir path
  Jobdir    = "#{Odir}/batch"
  # PrePkgdir = "#{Odir}/prefilter/refpkg"    # "#{Odir}/prefilter/refpkg/refpkg_info.tsv
  PreQuedir = "#{Odir}/prefilter/query"     # "#{Odir}/prefilter/query/#{que[:name]}/*.fa" (use .fa extension)
  PreFildir = "#{Odir}/prefilter/hmmsearch" # "#{Odir}/prefilter/hmmsearch/#{que[:name]}/{evalues.tsv,*.domtableout}"
  Cnkdir    = "#{Odir}/chunks"              # "#{Odir}/chunks/#{pkg[:name]}/{chunk,alignment,placement}"
  Resdir    = "#{Odir}/result"              # "#{Odir}/result/refpkg/#{pkg[:name]}/{all,each}/{query,alignment,placement,assign,extract}/#{que[:name]}"
  Logdir    = "#{Odir}/log/tasks"

  ## variables
  $fques    = []
  $rpkgs    = []
  $ex_lvs   = []

  ## validate evalue
  raise("#{Errmsg} -e, --evalue '#{ENV["evalue"]}' should be positive number.") if Evalue <= 0

  ## validate extract_levels
  ex_lvs_errmsg = "#{Errmsg} --extract-levels '#{Ex_lvs}' should be either -1, 0, positive integer, or comma separeted positive integers."
  case Ex_lvs
  when "-1", /^[0-9,]+$/
    a = Ex_lvs.split(",")
    raise(ex_lvs_errmsg) if a.size > 1 and a.include?("0")
    $ex_lvs = a.map(&:to_i).uniq.sort
  else raise(ex_lvs_errmsg)
  end

  ## Odir exist?
  $stderr.puts "\n\n#{Warmsg} output directory #{Odir} already exists. Overwrite it.\n\n" if OdirExist == "true"

  ### run
  NumStep  = Tasks.size
  Tasks.each.with_index(1){ |task, idx|
    Rake::Task[task].invoke(idx)
  }
end
# }}} default (run all tasks)


# {{{ tasks
desc "01-1a-A.validate_query"
task "01-1a-A.validate_query", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []
  script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  ## validate query and make a copy
  fques = Fques.split(",").inject([]){ |a, path| a += Dir[path.gsub("~", ENV["HOME"])].sort }
  fques.map{ |fque|
    if File.zero?(fque)
      $stderr.puts "#{Warmsg} #{fque} is empty. Skip it."
      next nil
    else
      File.open(fque){ |fr|
        if fr.gets[0] != ">"
          $stderr.puts "#{Warmsg} #{fque} is not valid fasta format. The first character should be '>'."
          next nil
        end
      } 
    end
    fque ## valid fasta
  }.compact

  $stderr.puts ["", "", "\e[1;32m===== check query file (N=#{fques.size}) \e[0m"]
  raise("#{Errmsg} no query file detected.") if fques.size == 0

  ### file name duplication check
  names = {}
  fques.each{ |fque|
    name = File.basename(fque).split(".")[0..-2]*"."
    raise("#{Errmsg} file name #{name} is not unique.") if names[name]
    names[name] = 1
  }

  mkdir_p PreQuedir
  idx = 0
  $fjsns = []
  fques.each{ |fque|
    ### parse fasta and check sequence length
    idx += 1
    name = File.basename(fque).split(".")[0..-2]*"."
    fa   = "#{PreQuedir}/#{name}.fa"
    fjsn = "#{fa}.json"

    if !File.exist?(fa) or !File.exist?(fjsn)
      outs << "ruby #{script} #{idx} #{MinSeqLen} #{name} #{fque} #{fa} #{fjsn}"
    end

    $fjsns << fjsn
  }

  next if outs.size == 0

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-1a-B.parse_query_info"
task "01-1a-B.parse_query_info", ["step"] do |t, args|
  $fques = []
  $fjsns.each{ |fjsn|
    sjsn = IO.readlines(fjsn)[0]
    $fques << eval(sjsn)
  }
end
desc "01-1b.validate_refpkg"
task "01-1b.validate_refpkg", ["step"] do |t, args|
  ## valdiate refpkg (do not make copy)
  ### [2022-06-15] sort files
  rpkgs = Rpkgs.split(",").sort_by{ |path| File.basename(path) }.inject([]){ |a, path| a += Dir[path.gsub("~", ENV["HOME"])].sort }
  $stderr.puts ["", "", "\e[1;32m===== check refpkg (N=#{rpkgs.size}) \e[0m"]
  raise("#{Errmsg} no refpkg directory detected.") if rpkgs.size == 0

  idx = 0
  names = {}
  $rpkgs = rpkgs.inject([]){ |a, rpkg|
    idx += 1
    fa    = Dir["#{rpkg}/*.fa"] + Dir["#{rpkg}/*.mfa"] + Dir["#{rpkg}/*.fasta"]
    fhmm  = Dir["#{rpkg}/*.hmm"]
    ftax  = Dir["#{rpkg}/taxon.tsv"] ;   ftax  = ftax.size == 0 ? nil : ftax[0]
    fpos  = Dir["#{rpkg}/position.tsv"]; fpos  = fpos.size == 0 ? nil : fpos[0]

    raise("#{Errmsg} #{rpkg} is not a directory.") unless File.directory?(rpkg)
    raise("#{Errmsg} #{rpkg} does not contain fasta file. #{rpkg}/*{.fa,.mfa,.fasta} should exist.") if fa.size == 0
    raise("#{Errmsg} #{rpkg} contains multiple fasta files. #{rpkg}/*{.fa,.mfa,.fasta} should be only one.") if fa.size > 1
    raise("#{Errmsg} #{rpkg} does not contain .hmm file. #{rpkg}/*.hmm should exist.") if fhmm.size == 0
    raise("#{Errmsg} #{rpkg} contains multiple .hmm files. #{rpkg}/*.hmm should be only one.") if fhmm.size > 1
    raise("#{Errmsg} #{fhmm[0]} contains .hmm file with multiple hmms. HMM should be only one per a file.") if IO.read(fhmm[0]).split(/\/\/\s+/).size > 1

    ### perse LENG of fhmm
    hmmlen = "0"
    open(fhmm[0]){ |fr|
      while l = fr.gets
        if l =~ /^LENG\s+(\d+)/
          hmmlen = $1
          break
        end
      end
    }

    ### [!!!] TODO: sequence (not name) in tree should be nonredundant
    ### [!!!] TODO: sequence name in tree should be same as sequence name in alignment
    ### [!!!] TODO: ftax format check
    ### [!!!] TODO: fpos format check

    name = File.basename(rpkg)
    raise("#{Errmsg} refpkg name #{name} is not unique.") if names[name]
    names[name] = 1

    h = { idx: idx, name: name, refpkg: rpkg, faln: fa[0], hmm: fhmm[0], ftax: ftax, fpos: fpos, hmmlen: hmmlen }
    $stderr.puts h.inspect
    a << h
  }

  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-2a.hmmsearch"
task "01-2a.hmmsearch", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []

  $fques.each{ |que|
    next if que[:numseq] == 0
    odir = "#{PreFildir}/#{que[:name]}"; mkdir_p odir
    $rpkgs.each{ |pkg|
      db   = pkg[:hmm]
      fa   = que[:fasta]
      fout = "#{odir}/#{pkg[:name]}.out"
      flog = "#{odir}/#{pkg[:name]}.log"
      outs << "hmmsearch --cpu 1 -Z #{Z} --notextw -o #{fout} #{db} #{fa} >#{flog} 2>&1"
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-2b.parse_hmmsearch"
task "01-2b.parse_hmmsearch", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []
  script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  pkgs    = $rpkgs.map{ |pkg| pkg[:name] }
  hmmlens = $rpkgs.map{ |pkg| pkg[:hmmlen] }

  $fques.each{ |que|
    fa    = que[:fasta]
    idir  = "#{PreFildir}/#{que[:name]}"
    next if Dir["#{idir}/*.out"].size == 0

    ## output:
    ## "#{Resdir}/refpkg/#{pkg[:name]}/each/query/{full_length,trimmed}/#{que[:name]}.fa" 
    ## "#{idir}/evalues.tsv" (best hit sequences for phylogenetic placement)
    outs << "ruby #{script} #{Trim_opt} #{Evalue} #{MinHmmLen} #{MinHmmLenFrc} #{fa} #{idir} #{Resdir} #{hmmlens*","} #{pkgs*" "}"
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-2c.merge_fasta"
task "01-2c.merge_fasta", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  outs   = []

  $rpkgs.each{ |pkg|
    odir = "#{Resdir}/refpkg/#{pkg[:name]}/all/query"; mkdir_p odir unless File.directory?(odir)

    ## merge trimmed position info
    fas = $fques.map{ |que| Dir["#{Resdir}/refpkg/#{pkg[:name]}/each/query/trimmed/#{que[:name]}.position"] }.flatten
    next if fas.size == 0
    fout = "#{odir}/trimmed.position"

    # outs << "cat #{fas*" "} >#{fout}" ## error if command line is too long
    open(fout, "w"){ |fw| fas.each{ |fa| fw.puts IO.read(fa) } }

    ## merge fasta
    %w|full_length trimmed|.each{ |type|
      fas = $fques.map{ |que| Dir["#{Resdir}/refpkg/#{pkg[:name]}/each/query/#{type}/#{que[:name]}.fa"] }.flatten
      next if fas.size == 0
      ## output:
      ## "#{Resdir}/refpkg/#{pkg[:name]}/each/all/{full_length,trimmed}.fa" 
      fout = "#{odir}/#{type}.fa"

      # outs << "cat #{fas*" "} >#{fout}" ## error if command line is too long
      open(fout, "w"){ |fw| fas.each{ |fa| fw.puts IO.read(fa) } }
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-3a.chunkify"
task "01-3a.chunkify", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa prepare chunkify --threads 1 --allow-file-overwriting"

  $rpkgs.each{ |pkg|
    # fas = "#{Resdir}/refpkg/#{pkg[:name]}/each/query/*/*.fa" ## hmmsearch alinged region
    fas = "#{Resdir}/refpkg/#{pkg[:name]}/each/query/trimmed/*.fa" ## hmmsearch alinged region
    next if Dir[fas].size == 0

    odir = "#{Cnkdir}/#{pkg[:name]}/chunk"; mkdir_p odir
    flog = "#{odir}/chunkify.log"
    outs << "#{cmd} --chunk-size #{C_size} --chunks-out-dir #{odir} --abundances-out-dir #{odir} --fasta-path #{fas} >#{flog} 2>&1"
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-3b.mafft_add"
task "01-3b.mafft_add", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

  script = "#{File.dirname(__FILE__)}/script/U2X.rb" ## "U" --> "X" in aligned fasta

  $rpkgs.each{ |pkg|
    fas = Dir["#{Cnkdir}/#{pkg[:name]}/chunk/chunk_*.fasta"].sort_by{ |i| File.basename(i).gsub(/^chunk_/, "").gsub(/\.fasta$/, "").to_i }
    next if fas.size == 0

    ### [!!!] faln should be nonredundant (sequences are same as tree) --> TODO: validation process
    faln = pkg[:faln]

    fas.each{ |fa|
      chnk  = File.basename(fa).gsub(/\.fasta$/, "") ## chunk_0, chunk_1, ...
      odir  = "#{Cnkdir}/#{pkg[:name]}/alignment/#{chnk}"; mkdir_p odir
      flog  = "#{odir}/mafft_add.log"
      fout  = "#{odir}/mafft_add.fa"
      ftmp  = "#{odir}/mafft_add.U2X.fa" ## used only by pplacer ('U' --> 'X')

      add   = "--add" ## --add or --addfragments

      ## [!!!] --maxiterate > 0 is not compatible with --keeplength
      if M_mafft == "FFT-NS-2" 
        option = "--retree 2"  ## fast and rough
        add    = "--addfragments"
      elsif M_mafft == "FFT-NS-i" 
        option = "--maxiterate 2"  ## fast and rough
        add    = "--addfragments"
      elsif M_mafft == "E-INS-i"
        option = "--genafpair" ## slow and accurate
        add    = "--add"
      else
        raise("#{Errmsg} --mafft-method #{M_mafft} is not available.")
      end

      out   = []
      out  << "mafft #{option} --anysymbol --thread 1 #{add} #{fa} --keeplength #{faln} >#{fout} 2>#{flog}" ## compatible with "U"
      out  << "ruby #{script} #{fout} #{ftmp}"
      outs << out*" && "

      ## [?] it could be automatic selection if option == ""
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-3c.pplacer"
task "01-3c.pplacer", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

  $rpkgs.each{ |pkg|
    fas = Dir["#{Cnkdir}/#{pkg[:name]}/alignment/chunk_*/mafft_add.U2X.fa"].sort_by{ |fa| fa.split("/")[-2].gsub(/^chunk_/, "").to_i }

    fas.each{ |fa|
      chnk  = fa.split("/")[-2]
      fpkg  = pkg[:refpkg] ### [!!!] should be nonredundant (sequences are same as tree)
      odir  = "#{Cnkdir}/#{pkg[:name]}/placement/#{chnk}"; mkdir_p odir
      flog  = "#{odir}/pplacer.log"

      outs << "pplacer -j 1 --out-dir #{odir} -c #{fpkg} #{fa} >#{flog} 2>&1"
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-3d.unchunkify"
task "01-3d.unchunkify", ["step"] do |t, args|
  ### [!!!] [2020-11-17] Use of parallel cause "core dump". DO NOT USE parallel
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa prepare unchunkify --threads 1 --allow-file-overwriting"

  script = "#{File.dirname(__FILE__)}/script/merge_jplace.rb"

  $rpkgs.each{ |pkg|
    fplcs = "#{Cnkdir}/#{pkg[:name]}/placement/chunk_*/*.jplace" ## mafft_add.jplace
    next if Dir[fplcs].size == 0

    fabus = "#{Cnkdir}/#{pkg[:name]}/chunk/abundances_*.json"
    odir  = "#{Resdir}/refpkg/#{pkg[:name]}/each/placement"; mkdir_p odir ## #{Resdir}/refpkg/#{pkg[:name]}/each/placement/*.jplace
    flog  = "#{odir}/unchunkify.log"

    out   = []
    out  << "#{cmd} --out-dir #{odir} --abundances-path #{fabus} --jplace-path #{fplcs} >#{flog} 2>&1"

    ### merge the all jplace file
    odir  = "#{Resdir}/refpkg/#{pkg[:name]}/all/placement"; mkdir_p odir
    fall  = "#{odir}/all.jplace"

    out  << "ruby #{script} #{Resdir}/refpkg/#{pkg[:name]}/each/placement/*.jplace >#{fall}"
    outs << out*" && "
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)

  ### [!!!] [2020-11-17] DO NOT USE PARALLEL
  # RunBatch.call(t, Jobdir, Ncpu, Logdir)
  RunBatch.call(t, Jobdir, 1, Logdir)
end
desc "01-3e.unchunkify_alignment"
task "01-3e.unchunkify_alignment", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

  script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  $rpkgs.each{ |pkg|
    fas  = "#{Cnkdir}/#{pkg[:name]}/chunk/chunk_*.fasta" ## trimmed seq, hashed 
    alns = "#{Cnkdir}/#{pkg[:name]}/alignment/chunk_*/mafft_add.fa" ## aligned seq, hashed
    ques = $fques.map{ |que| Dir["#{Resdir}/refpkg/#{pkg[:name]}/each/query/trimmed/#{que[:name]}.fa"] }.flatten ## trimmed seq, labeled 
    next if ques.size == 0

    ### write unchunked alignment in #{Resdir}/refpkg/#{pkg[:name]}/all/alignment/aligned.fa
    odir  = "#{Resdir}/refpkg/#{pkg[:name]}/all/alignment"; mkdir_p odir unless File.directory?(odir)
    faall = "#{odir}/aligned.fa"

    ### write unchunked alignment in #{Resdir}/refpkg/#{pkg[:name]}/each/alignment/#{que[:name]}/aligned.fa
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}/each/alignment"; mkdir_p bdir

    outs << "ruby #{script} #{pkg[:faln]} '#{fas}' '#{alns}' '#{ques*","}' #{faall} #{bdir}"
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4a.info"
task "01-4a.info", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa examine info --threads 1 --allow-file-overwriting"

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all (merged) query files
    fplcs = "#{bdir}/all/placement/all.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size == 0

    odir  = "#{bdir}/all/placement"; mkdir_p odir
    flog  = "#{odir}/all.info"
    outs << "#{cmd} --jplace-path #{fplcs} >#{flog} 2>&1"

    ### each query file
    $fques.each{ |que|
      fplcs = "#{bdir}/each/placement/#{que[:name]}.jplace" ## "#{que[:name]}.jplace"
      next unless File.exist?(fplcs)

      odir  = "#{bdir}/each/placement"
      flog  = "#{odir}/#{que[:name]}.info"
      outs << "#{cmd} --jplace-path #{fplcs} >#{flog} 2>&1"
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4b.lwr"
task "01-4b.lwr", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa examine lwr --threads 1 --allow-file-overwriting"

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all (merged) query files
    fplcs = "#{bdir}/all/placement/all.jplace"
    # fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size == 0

    odir  = "#{bdir}/all/placement"; mkdir_p odir
    flog  = "#{odir}/all.lwr.log"
    pref  = "all.lwr_"
    outs << "#{cmd} --out-dir #{odir} --file-prefix #{pref} --jplace-path #{fplcs} >#{flog} 2>&1"

    ### each query file
    $fques.each{ |que|
      fplcs = "#{bdir}/each/placement/#{que[:name]}.jplace" ## "#{que[:name]}.jplace"
      next unless File.exist?(fplcs)

      odir  = "#{bdir}/each/placement"
      flog  = "#{odir}/#{que[:name]}.lwr.log"
      pref  = "#{que[:name]}.lwr_"
      outs << "#{cmd} --out-dir #{odir} --file-prefix #{pref} --jplace-path #{fplcs} >#{flog} 2>&1"
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4c.edpl"
task "01-4c.edpl", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa examine edpl --no-list-file --threads 1 --allow-file-overwriting"
  ### [!!!] --no-list-file: need to limit memory usage.
  ### If set, do not write out the EDPL per pquery, but just the histogram file. As the list needs to keep all pquery names in memory (to get the correct order), the memory requirements might be too large. In that case, this option can help.

  ## [!!!] memory requirement is very high even when --no-list-file is used.
  $stderr.puts "Currently, do not execute 'gappa examine edpl' due to high memory requirement."
  next ## SKIP THIS TASK !!!

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all query files
    fplcs = "#{bdir}/all/placement/all.jplace"
    # fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size == 0

    odir  = "#{bdir}/all/placement"; mkdir_p odir
    flog  = "#{odir}/all.edpl.log"
    pref  = "all.edpl_"
    outs << "#{cmd} --out-dir #{odir} --file-prefix #{pref} --jplace-path #{fplcs} >#{flog} 2>&1"

    ### each query file
    $fques.each{ |que|
      fplcs = "#{bdir}/each/placement/#{que[:name]}.jplace" ## "#{que[:name]}.jplace"
      next unless File.exist?(fplcs)

      odir  = "#{bdir}/each/placement"
      flog  = "#{odir}/#{que[:name]}.edpl.log"
      pref  = "#{que[:name]}.edpl_"
      outs << "#{cmd} --out-dir #{odir} --file-prefix #{pref} --jplace-path #{fplcs} >#{flog} 2>&1"
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4d.assign"
task "01-4d.assign", ["step"] do |t, args|
  ### option might be changed between v0.6.0 and v0.7.1
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  # cmd  = "gappa examine assign --threads 1 --per-query-results --krona --allow-file-overwriting" ### for v0.7.1
  cmd  = "gappa examine assign --threads 1 --krona --allow-file-overwriting" ### for v0.6.0

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"
    ftax = pkg[:ftax]

    ($stderr.puts "For #{pkg}, taxon.tsv is not found. Skip the taxonomy/clade assignment step."; next) unless ftax

    ### all query files
    fplcs = "#{bdir}/all/placement/all.jplace"
    # fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size == 0

    odir  = "#{bdir}/all/assign"; mkdir_p odir
    flog  = "#{odir}/assign.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} --taxon-file #{ftax} >#{flog} 2>&1"

    ### each query file
    $fques.each{ |que|
      fplcs = "#{bdir}/each/placement/#{que[:name]}.jplace" ## "#{que[:name]}.jplace"
      next unless File.exist?(fplcs)

      odir  = "#{bdir}/each/assign/#{que[:name]}"; mkdir_p odir
      flog  = "#{odir}/assign.log"
      outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} --taxon-file #{ftax} >#{flog} 2>&1"
    }
  }

  # next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4e.extract"
task "01-4e.extract", ["step"] do |t, args|
  next if $ex_lvs == [-1] ### do not extract

  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa prepare extract --threads 1 --allow-file-overwriting"

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"
    ftax = pkg[:ftax]

    ($stderr.puts "For #{pkg}, taxon.tsv is not found. Skip the taxonomy/clade assignment step."; next) unless ftax

    ### parse levels
    max_lv = IO.readlines(ftax).map{ |l| i = l.chomp.split("\t")[1]; i ? i.split(/\s*;\s*/).size : 0 }.sort[-1]
    next if max_lv == 0
    lvs = if $ex_lvs == [0] then (1..max_lv).to_a
    else $ex_lvs.select{ |i| i <= max_lv }
    end

    lvs.each{ |lv|
      ### all query files
      fplcs = "#{bdir}/all/placement/all.jplace"
      # fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
      next if Dir[fplcs].size == 0

      fas   = "#{bdir}/all/query/full_length.fa"
      odir  = "#{bdir}/all/extract/level#{lv}"; mkdir_p odir
      flog  = "#{odir}/extract.log"
      fsvg  = "#{odir}/color_tree.svg"

      fcld  = "#{odir}/clade.tsv"
      open(fcld, "w"){ |fw|
        fw.puts IO.readlines(ftax).map{ |l| a = l.chomp.split("\t"); [a[0], a[1].split(/\s*;\s*/)[lv - 1]]*"\t" }
      }

      outs << "#{cmd} --samples-out-dir #{odir} --sequences-out-dir #{odir} --jplace-path #{fplcs} --fasta-path #{fas} --clade-list-file #{fcld} >#{flog} 2>&1"

      ### each query file
      $fques.each{ |que|
        fplcs = "#{bdir}/each/placement/#{que[:name]}.jplace" ## "#{que[:name]}.jplace"
        next unless File.exist?(fplcs)

        fas   = "#{bdir}/each/query/full_length/#{que[:name]}.fa"
        odir  = "#{bdir}/each/extract/#{que[:name]}"; mkdir_p odir
        flog  = "#{odir}/extract.log"
        fsvg  = "#{odir}/color_tree.svg"
        fcld  = "#{bdir}/all/extract/level#{lv}/clade.tsv" ## written above
        outs << "#{cmd} --samples-out-dir #{odir} --sequences-out-dir #{odir} --jplace-path #{fplcs} --fasta-path #{fas} --clade-list-file #{fcld} >#{flog} 2>&1"
      }
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4f.graft"
task "01-4f.graft", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa examine graft --threads 1 --fully-resolve --name-prefix Q_ --allow-file-overwriting"

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all query files
    fplcs = "#{bdir}/all/placement/all.jplace"
    # fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size == 0

    odir  = "#{bdir}/all/graft"; mkdir_p odir
    flog  = "#{odir}/graft.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"

    ### each query file
    $fques.each{ |que|
      fplcs = "#{bdir}/each/placement/#{que[:name]}.jplace" ## "#{que[:name]}.jplace"
      next unless File.exist?(fplcs)

      odir  = "#{bdir}/each/graft/#{que[:name]}"; mkdir_p odir
      flog  = "#{odir}/graft.log"
      outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4g.heat-tree"
task "01-4g.heat-tree", ["step"] do |t, args|
  ### option might be changed between v0.6.0 and v0.7.1
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  opt_svg = "--svg-tree-shape circular --color-list viridis --reverse-color-list --svg-tree-stroke-width 3 --svg-tree-ladderize"
  prefix  = "tree"
  cmd  = "gappa examine heat-tree --threads 1 #{opt_svg} --tree-file-prefix #{prefix} --write-newick-tree --write-nexus-tree --write-phyloxml-tree --write-svg-tree --allow-file-overwriting"

  script = "#{File.dirname(__FILE__)}/script/nexus2itol.rb"

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all query files
    fplcs = "#{bdir}/all/placement/all.jplace"
    # fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size == 0

    odir  = "#{bdir}/all/heat-tree"; mkdir_p odir
    flog  = "#{odir}/heat-tree.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1 && ruby #{script} #{odir}/#{prefix}.nexus #{odir}/#{prefix}.itol.txt"

    ### each query file
    $fques.each{ |que|
      fplcs = "#{bdir}/each/placement/#{que[:name]}.jplace" ## "#{que[:name]}.jplace"
      next unless File.exist?(fplcs)

      odir  = "#{bdir}/each/heat-tree/#{que[:name]}"; mkdir_p odir
      flog  = "#{odir}/heat-tree.log"
      outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1 && ruby #{script} #{odir}/#{prefix}.nexus #{odir}/#{prefix}.itol.txt"
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4h.aligned_position"
task "01-4h.aligned_position", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

  script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  $rpkgs.each{ |pkg|
    bdir    = "#{Resdir}/refpkg/#{pkg[:name]}"
    fpos    = pkg[:fpos]
    frefaln = pkg[:faln]

    ($stderr.puts "For #{pkg}, position.tsv is not found. Skip to check alignment positions."; next) unless fpos

    ### all query files
    fplcs = "#{bdir}/all/placement/all.jplace"
    # fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size == 0

    odir    = "#{bdir}/all/alignment"
    fallaln = "#{odir}/aligned.fa" ## should exist.
    fasn    = "#{bdir}/all/assign/per_query.tsv" ## use result of assignment. If not exist, do not use assignment result
    flog    = "#{odir}/aligned_position.log"
    outs << "ruby #{script} #{odir} #{fpos} #{fallaln} #{frefaln} #{fasn} >#{flog} 2>&1"

    ### each query file
    $fques.each{ |que|
      fplcs = "#{bdir}/each/placement/#{que[:name]}.jplace" ## "#{que[:name]}.jplace"
      next unless File.exist?(fplcs)

      odir    = "#{bdir}/each/alignment/#{que[:name]}"
      fallaln = "#{odir}/aligned.fa"                ## should exist.
      fasn    = "#{bdir}/each/assign/#{que[:name]}/per_query.tsv" ## use result of assignment. If not exist, do not use assignment result
      flog    = "#{odir}/aligned_position.log"
      outs << "ruby #{script} #{odir} #{fpos} #{fallaln} #{frefaln} #{fasn} >#{flog} 2>&1"
    }
  }

  # next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-5a.krd"
task "01-5a.krd", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa analyze krd --threads 1 --allow-file-overwriting"

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all query files
    fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size < 2

    odir  = "#{bdir}/all/krd"; mkdir_p odir
    flog  = "#{odir}/krd.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-5b.edgepca"
task "01-5b.edgepca", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa analyze edgepca --threads 1 --write-nexus-tree --write-phyloxml-tree --write-svg-tree --allow-file-overwriting"
  ## newick tree is same as tree in refpkg (no color information) --> do not use --write-newick-tree

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all query files
    fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size < 2

    odir  = "#{bdir}/all/edgepca"; mkdir_p odir
    flog  = "#{odir}/edgepca.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-5c.squash"
task "01-5c.squash", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa analyze squash --threads 1 --write-nexus-tree --write-phyloxml-tree --write-svg-tree --allow-file-overwriting"
  ## newick tree is same as tree in refpkg (no color information) --> do not use --write-newick-tree
  ## [!!!] this task will generate (2n - 2) * 3 tree files where n is num of input query fasta files.

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all query files
    fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size < 2

    odir  = "#{bdir}/all/squash"; mkdir_p odir
    flog  = "#{odir}/squash.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-5d.dispersion"
task "01-5d.dispersion", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []
  cmd  = "gappa analyze dispersion --threads 1 --write-nexus-tree --write-phyloxml-tree --write-svg-tree --allow-file-overwriting"
  ## newick tree is same as tree in refpkg (no color information) --> do not use --write-newick-tree

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all query files
    fplcs = "#{bdir}/each/placement/*.jplace" ## "#{que[:name]}.jplace"
    next if Dir[fplcs].size < 2

    odir  = "#{bdir}/all/dispersion"; mkdir_p odir
    flog  = "#{odir}/dispersion.log"
    outs << "#{cmd} --out-dir #{odir} --mass-norm absolute --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-6a.aa_feature"
task "01-6a.aa_feature", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

  script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  $rpkgs.each{ |pkg|
    bdir = "#{Resdir}/refpkg/#{pkg[:name]}"

    ### all query files
    fa    = "#{bdir}/all/query/trimmed.fa"
    next unless File.exist?(fa)

    odir  = "#{bdir}/all/feature/aa"; mkdir_p odir
    fout  = "#{odir}/aa_feature.tsv"
    flog  = "#{odir}/aa_feature.log"
    outs << "ruby #{script} #{fa} #{fout} >#{flog} 2>&1"

    ### each query file
    $fques.each{ |que|
      fa    = "#{bdir}/each/query/trimmed/#{que[:name]}.fa"
      next unless File.exist?(fa)

      odir  = "#{bdir}/each/feature/aa"; mkdir_p odir unless File.directory?(odir)
      fout  = "#{odir}/aa_feature.tsv"
      flog  = "#{odir}/aa_feature.log"
      outs << "ruby #{script} #{fa} #{fout} >#{flog} 2>&1"
    }
  }

  next if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit")
  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}} tasks
