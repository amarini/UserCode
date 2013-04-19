

.PHONY: BaseAnalysis
BaseAnalysis:
	cd BaseAnalysis && make

.PHONY: ZAnalysis
ZAnalysis: BaseAnalysis
	cd ZAnalysis && make

.PHONY: GAnalysis
GAnalysis: BaseAnalysis
	cd GAnalysis && make

.PHONY: zexe
zexe:
	cd ZAnalysis && make zexe
	[ -h zexe ] || ln -s ZAnalysis/zexe

