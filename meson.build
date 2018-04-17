project('gc', ['c', 'cpp'],
	default_options : ['c_std=c11', 'cpp_std=c++1z'])
compiler = meson.get_compiler('cpp')

basic_lib = static_library('basic', ['src/intstack.cpp',
				     'src/graph.cpp',
				     'src/options.cpp',
				     'src/statistics.cpp'])
basic_dep = declare_dependency(link_with:basic_lib,
			       include_directories: include_directories('src'))

minicsp_inc = include_directories('minicsp')
minicsp_lib = library('minicsp', ['minicsp/core/solver.cpp',
				  'minicsp/core/cons.cpp',
				  'minicsp/core/utils.cpp'])
minicsp_dep = declare_dependency(link_with: minicsp_lib,
				 include_directories: minicsp_inc)

executable('gc', ['src/gc.cpp',
		  'src/prop.cpp',
		  'src/rewriter.cpp',
		  'src/mycielski.cpp'],
	   dependencies: [minicsp_dep, basic_dep])

use_gecode = get_option('use_gecode')
if use_gecode
  gecode_root = get_option('gecode')
  gecode_libs = gecode_root+'/lib'
  gecode_inc = include_directories(gecode_root+'/include')
  gecodedriver_lib = compiler.find_library('gecodedriver', dirs:gecode_root+'/lib')
  gecodesearch_lib = compiler.find_library('gecodesearch', dirs:gecode_root+'/lib')
  gecodeminimodel_lib = compiler.find_library('gecodeminimodel', dirs:gecode_root+'/lib')
  gecodeint_lib = compiler.find_library('gecodeint', dirs:gecode_root+'/lib')
  gecodeset_lib = compiler.find_library('gecodeset', dirs:gecode_root+'/lib')
  gecodefloat_lib = compiler.find_library('gecodefloat', dirs:gecode_root+'/lib')
  gecodekernel_lib = compiler.find_library('gecodekernel', dirs:gecode_root+'/lib')
  gecodesupport_lib = compiler.find_library('gecodesupport', dirs:gecode_root+'/lib')
  
  gecode_dep = declare_dependency(dependencies: [gecodedriver_lib,
						 gecodesearch_lib,
						 gecodeminimodel_lib,
						 gecodeint_lib,
						 gecodeset_lib,
						 gecodefloat_lib,
						 gecodekernel_lib,
						 gecodesupport_lib],
				  include_directories: gecode_inc)
  
  executable('gecode-gc', ['sota/graph-color.cpp'],
	     dependencies: [gecode_dep, basic_dep])
endif

executable('minicolor-converter',
	   ['sota/converter/main.c',
	    'sota/converter/converter.c'])

gzip_lib = compiler.find_library('z')
gzip_dep = declare_dependency(dependencies: gzip_lib)
executable('minicolor',
	   ['sota/MiniMerge/core/Main.cpp',
	    'sota/MiniMerge/core/Solver.cpp'],
	   include_directories: include_directories(['sota/MiniMerge/mtl']),
	   dependencies: gzip_dep)

