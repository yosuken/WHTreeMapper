
require 'rake'

Trim_opt, e_thre, minhmmlen, fa, idir, resdir = ARGV
E_thre = e_thre.to_f
MinHmmLen = minhmmlen.to_i

# {{{ range function
class Range
  include Comparable

  def <=>(other)
    self.min <=> other.min
  end
  def overlap?(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return false
    else
      return true
    end
  end
  def &(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return nil
    else
      return ([s.min, t.min].max)..([s.max, t.max].min)
    end
  end
  def |(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return [s, t]
    else
      return [([s.min, t.min].min)..([s.max, t.max].max)]
    end
  end
  def include?(t)
    s = self
    o = s & t  ## overlap of s and t
    return false unless o
    if o.min == t.min and o.max == t.max ## s includes t
      return true
    else
      return false
    end
  end
end
def merge_ranges(ranges) ## ranges = [3..300, 320..500, 504..732, ...]
  return ranges if ranges.size < 2

  rs  = ranges.sort
  (-rs.size..-2).each{ |j| ## index from right
    merged = rs[j] | rs[j+1]
    case merged.size
    when 1 ## overlap detected
      rs[j] = merged[0]
      rs.delete_at(j+1)
    when 2 ## overlap not detected
      ## do nothing
    else raise
    end
  }
  return rs
end
# }}} range function

fins = Dir["#{idir}/*.out"]
qname = File.basename(fa).split(".")[0..-2]*"." ## query name

evalues    = Hash.new{ |h, i| h[i] = Array.new(fins.size, "-") } ## store full i-Evalue
domtblout  = Hash.new{ |h, i| h[i] = Array.new(fins.size, nil) }
hmm_coords = Hash.new{ |h, i| h[i] = Array.new(fins.size, nil) } ## hit hmm coordinates
ali_coords = Hash.new{ |h, i| h[i] = Array.new(fins.size, nil) } ## hit ali coordinates
labs       = []

fins.each.with_index{ |fin, idx|
  labs << File.basename(fin).gsub(/\.out$/, "") ## refpkg label
  n_ent = 0
  flag  = ""
  gid   = ""

# Query:       1-cysPrx_C  [M=40]
# Accession:   PF10417.9
# Description: C-terminal domain of 1-Cys peroxiredoxin
# Scores for complete sequences (score includes all domains):
#    --- full sequence ---   --- best 1 domain ---    -#dom-
#     E-value  score  bias    E-value  score  bias    exp  N  Sequence  Description
#     ------- ------ -----    ------- ------ -----   ---- --  --------  -----------
#     4.6e-12   51.2   0.4    6.3e-12   50.7   0.4    1.2  1  gene1
#       0.004   22.6   0.5      0.016   20.6   0.0    2.4  1  geneX
#   ------ inclusion threshold ------
#       0.011   21.2   0.0       0.03   19.7   0.0    1.8  1  geneY

# >> ERS488460_ERR1719473_1814004
#    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
#  ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#    1 ?    6.6   0.0    0.0083   2.1e+02      90     135 ..      16      57 ..       5      74 .. 0.79
#    2 ?    7.7   0.0     0.004     1e+02      89     134 ..     123     164 ..      94     173 .. 0.80
#
  open(fin){ |fr|
    while l = fr.gets
      if l =~ /^Query:/  ## detect new entry
        n_ent += 1
        raise("multiple entries detected. abort.") if n_ent > 1
        flag = "parse_full"
      elsif l =~ /^$/ ## reset parse flag if a blank line comes
        flag = ""
        gid  = ""
      elsif n_ent == 1 and flag == "parse_full" and l =~ /^\s+\d/
        evalue, gid = l.strip.split(/\s+/).values_at(0, 8) ## full i-evalue
        next if evalue.to_f >= E_thre
        evalues[gid][idx] = evalue.to_f
      elsif l =~ /^>> (\S+)/ ## reset parse flag if a blank line comes
        gid  = $1
        next if evalues[gid][idx] == "-" ## full i-evalue is above threshold ==> do not parse
        flag = "parse_each_domain"
      elsif n_ent == 1 and flag == "parse_each_domain" and l =~ /^\s+\d/
        evalue, hmm_fm, hmm_to, ali_fm, ali_to = l.strip.split(/\s+/).values_at(5, 6, 7, 9, 10) ## domain i-evalue
        next if evalue.to_f >= E_thre
        hmm_coords[gid][idx] ||= []
        hmm_coords[gid][idx] <<  (hmm_fm.to_i..hmm_to.to_i)
        ali_coords[gid][idx] ||= []
        ali_coords[gid][idx] <<  (ali_fm.to_i..ali_to.to_i)
      end
    end
  }
}

### check hmm hit length
evalues.each_key{ |gid|
  ary = evalues[gid]
  (0..ary.size-1).each{ |idx|
    evalue = ary[idx]
    next if evalue == "-"
    hmm_co = hmm_coords[gid][idx]
    ali_co = ali_coords[gid][idx]
    if hmm_co.nil?
      evalues[gid][idx] = "-"
    else
      if Trim_opt == "merge"
        hmm_co = merge_ranges(hmm_co)
        ali_co = merge_ranges(ali_co)
      elsif Trim_opt == "largest"
        lens = hmm_co.map{ |r| r.max - r.min + 1 }
        _idx = lens.each_with_index.max[1] ## take index of max value
        hmm_co = [hmm_co[_idx]] ## take the largest
        ali_co = [ali_co[_idx]] ## take the corresponding hit
      else raise("unknown trim_opt: #{Trim_opt}")
      end
      hmmlen = hmm_co.map{ |r| r.max - r.min + 1 }.inject(&:+)
      if hmmlen < MinHmmLen ## hmm hit length is not enough
        evalues[gid][idx]    = "-"
        hmm_coords[gid][idx] = nil
        ali_coords[gid][idx] = nil
      else
        hmm_coords[gid][idx] = hmm_co.map{ |r| "#{r.min}-#{r.max}" }*","
        ali_coords[gid][idx] = [ali_co[0].min, ali_co[-1].max]
      end
    end
  }
}

### evalues.tsv
lab2gids = Hash.new{ |h, i| h[i] = [] }
gid2ent  = {}
gid2pos  = {} ## 1-based position [from, to]
open("#{idir}/evalues.tsv", "w"){ |fw|
  fw.puts ["seq", labs]*"\t"

  IO.read(fa).split(/^>/)[1..-1].each{ |ent|
    _lab, *seq = ent.split("\n")
    gid = _lab.split(/\s+/)[0]

    a = evalues[gid] ## e.g., ["-", 3e-10, 5e-9, "-"]
    fw.puts [gid, a]*"\t"

    next if a.uniq == ["-"]
    a = a.map{ |i| i == "-" ? 100000 : i }
    a = a.map.with_index{ |i, idx| [i, idx] }
    minidx = a.sort_by{ |i, idx| i }[0][1]
    lab = labs[minidx]

    lab2gids[lab] << gid
    gid2ent[gid]   = ent

    ali_co = ali_coords[gid][minidx]
    raise unless ali_co
    gid2pos[gid] = ali_co
  }
}

lab2gids.each{ |lab, gids|
  # next if gids.size == 0

  odir1 = "#{resdir}/refpkg/#{lab}/each/query/full_length"; mkdir_p odir1 unless File.directory?(odir1)
  odir2 = "#{resdir}/refpkg/#{lab}/each/query/trimmed";     mkdir_p odir2 unless File.directory?(odir2)

  open("#{odir1}/#{qname}.fa", "w"){ |fw| gids.each{ |gid| fw.puts ">"+gid2ent[gid] } }
  fwp = open("#{odir2}/#{qname}.position", "w")
  open("#{odir2}/#{qname}.fa", "w"){ |fw| gids.each{ |gid|
    ent = gid2ent[gid]
    lab, *seq = ent.split("\n")
    seq = seq.join.gsub(/\s+/, "")
    pos = gid2pos[gid]
    len = pos[1].to_i - pos[0].to_i + 1
    raise if len > seq.size
    trimmed = seq[pos[0].to_i - 1, len]
    lab2 = [lab.strip, pos*"-"]*"\t"

    fw.puts [">"+lab, trimmed]
    fwp.puts lab2
  } }
  fwp.close
}

