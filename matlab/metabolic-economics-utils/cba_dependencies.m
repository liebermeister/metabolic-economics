function cba_dependencies()

if ~exist('mnt_version','file'),
  error('Please install the Metabolic Network Toolbox (https://github.com/wolframliebermeister/mnt)');
end

mnt_dependencies
