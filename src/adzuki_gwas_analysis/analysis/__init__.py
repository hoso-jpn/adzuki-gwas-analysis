"""Analysis functions for one already-validated GWAS summary-statistics dataset.

Everything under this package operates on a single dataset at a time -- never
more than one of the 6 Dryad files is held in memory simultaneously -- and
assumes its input has already passed
:func:`adzuki_gwas_analysis.validate.validate_dataset`. No function here
re-implements or loosens that validation; a caller that skips it gets whatever
:mod:`pandas` does with malformed data, which is why every CLI/pipeline entry
point in :mod:`adzuki_gwas_analysis.analysis.pipeline` runs validation first
and refuses to produce output if it fails.
"""

from __future__ import annotations
