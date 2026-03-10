
STDOUT.sync = true; STDERR.sync = true

odir, fpos, fallaln, frefaln, fasn = ARGV

## parse fpos
lab2pos = {}
IO.readlines(fpos).each.with_index(1){ |l, idx|
  next if l =~ /^#/ ## skip comment lines (begins with "#")

  lab, pos = l.chomp.split("\t")
  lab2pos[lab] = pos.split(",").map{ |i|
    if i.strip !~ /^[0-9]+$/
      raise("\e[1;31mError:\e[0m invalid position value in line #{idx} of #{fpos} - '#{l}'. postion should be positive integer (and comma separeted if multiple values)")
    end
    i.strip.to_i
  }
}

## parse fasn --- assign result (might not exist --> [!!!] in this case, reg2info is empty) 
reg2info = {} ## (name LWR fract aLWR afract taxopath) --> parse fract, taxpath
reg2set  = {} ## consider redundant sequence (ERR315859.153160.1_2_5_1 has redundant mapping result 'ERR315859.153160.1_2_5_1 ;ERR315859.32048204.1_3_6_1')
if File.exist?(fasn)
  IO.readlines(fasn)[1..-1].each{ |l|
    regs, fract, tax = l.chomp.split("\t").values_at(0, 2, 5)

    # regs = regs.split(";").map{ |reg| ## parse ';' separated ids
    #   reg.strip.split(/\s+/)[0] ## exclude additional info
    # }
    regs = regs.split(/\s*;/).map{ |reg| reg.strip }

    regs.each{ |reg|
      reg2info[reg] ||= [] ## store [fract, tax]
      reg2info[reg] << [fract, tax]

      reg2set[reg] = regs
    }
  }

  reg2info.each_key{ |reg|
    a = reg2info[reg]
    best = a.sort_by{ |i| - i[0].to_f }[0] ## extract if fract (fract of LWR) is higher
    reg2info[reg] = best ## only highest likelihood hit
  }
else 
  $stderr.puts "assign result: #{fasn} does not exist. generate result without taxon/clade assignment info."
end

## parse last reg in frefaln
ref_last = ""
IO.read(frefaln).split(/^>/)[1..-1].each{ |ent|
  lab, *_ = ent.split("\n")
  reg = lab.split(/\s+/)[0]
  ref_last = reg ## lastly, the last reg will be stored
}

## parse fallaln
open("#{odir}/aligned_position.tsv", "w"){ |fw|
  fw.puts ["query", lab2pos.keys, "fract", "taxpath"]*"\t"

  flag = 0 ## ref or que

  IO.read(fallaln).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split("\n")
    reg = lab.split(/\s+/)[0]
    seq = seq.join.gsub(/\s+/, "")

    if flag == 1
      out = [reg]

      ## parse amino acids in specified position
      lab2pos.each{ |lab, pos| out << pos.map{ |i| seq[i-1] }*"," }

      ## add taxon/clade info
      out << (reg2info[reg] ? reg2info[reg] : ["NA", "NA"])

      reg, *info = out
      fw.puts [reg, info]*"\t"

      # raise("reg: #{reg} - not found in reg2set[reg].") unless reg2set[reg]
      # reg2set[reg].each{ |_reg| ### make multiple lines if member of a set is multiple.
      #   fw.puts [_reg, info]*"\t"
      # }
    end

    flag = 1 if ref_last == reg ## make flag of ref --> que
  }
}
