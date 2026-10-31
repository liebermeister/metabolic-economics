function cba_dependencies()

if ~exist('mnt_version','file'),
  error('Please install the Metabolic Network Toolbox (https://github.com/liebermeister/metabolic-network-toolbox)');
end

if exist('package_version_check','file')
  package_version_check('The Metabolic Network Toolbox', 'mnt_version', '1.3.0');
end

mnt_dependencies
