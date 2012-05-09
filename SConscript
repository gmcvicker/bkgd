
Import('env', 'libgenome_db', 'libgenome', 'libanalysis')


env = env.Clone()

# add CFLAGS and LFLAGS for GSL library



# env.PrependENVPath(self, name, os.environ['PATH']})

env.ParseConfig(['gsl-config',  '--cflags', '--libs'])

env.MergeFlags({'LIBS' : ['analysis']})

objs = env.Object(['bkgd.c', 'bkgd_param.c', 'bkgd_interp.c',
                   'interp_tab.c', 'bkgd_intg.c',
                   'bkgd_point.c', 'bkgd_evo_mdl.c',
                   'bkgd_evo_mdl_param.c',
                   'cltype.c', 'branch.c', 'model.c', 'bkgd_data.c'])

cb_env = env.Clone()
p = cb_env.Program('calc_bkgd', ['calc_bkgd.c', objs])


p += env.Program('bkgd_mle', ['bkgd_mle.c', objs])

p += env.Program('bkgd_expect_sites', ['bkgd_expect_sites.c', objs])

p += env.Program('bkgd_subst_prob_report', ['bkgd_subst_prob_report.c', objs])

p += env.Program('bkgd_hist', ['bkgd_hist.c', objs])


Depends(p, [libgenome_db, libgenome, libanalysis])

