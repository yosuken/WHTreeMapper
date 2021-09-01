
### merge jplace files
### usage: ruby <this scrript> <jplace files>

### [!!! CAUTION]
### This script has strict assumption on input files (e.g.,  all jplace files has the same "fields" field). See below.
### 1. tree field should be identical as string.
### 2. fields field shuold be identical.

fins = ARGV

require 'json'

h = {}
fins.each{ |fin|
  ($stderr.puts "input file (#{fin}) is empty. skip it." ; next) if File.zero?(fin)

  str = IO.read(fin)
  jsn = JSON.parse(str, max_nesting: false)
  # $stderr.puts "parsed: #{fin}"

  if h.size == 0 ### first .jplace file
    h = jsn
  else
    ### Is tree identical?
    raise("tree of #{fin} is different from others. Aborting.") if  h["tree"] != jsn["tree"]
    ### Is fields identical?
    raise("fields of #{fin} is different from others. Even the order shuold be identical. Aborting.") if  h["fields"] != jsn["fields"]
    ### add placements
    h["placements"] += jsn["placements"]
  end
}

puts JSON.pretty_generate(h)
# puts JSON.generate(jsn)
