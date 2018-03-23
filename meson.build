project('gc', 'cpp',
	default_options : ['c_std=c11', 'cpp_std=c++1z'])

minicsp_inc = include_directories('minicsp')
minicsp_lib = library('minicsp', ['minicsp/core/solver.cpp',
				  'minicsp/core/cons.cpp',
				  'minicsp/core/utils.cpp'])
minicsp_dep = declare_dependency(link_with: minicsp_lib,
				 include_directories: minicsp_inc)
  
executable('gc', ['src/gc.cpp',
		  'src/intstack.cpp',
		  'src/graph.cpp',
		  'src/prop.cpp',
		  'src/options.cpp',
			'src/statistics.cpp'],
	   dependencies: minicsp_dep)
