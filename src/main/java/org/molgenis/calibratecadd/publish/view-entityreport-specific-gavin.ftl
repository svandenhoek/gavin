<#-- modal header, view-entityreport-specific-gavin_ruleguide.ftl -->
<div class="modal-header">
    <button type="button" class="close" data-dismiss="modal" aria-hidden="true">&times;</button>
    <h4 class="modal-title">DataSet: ${entityMetadata.getLabel()?html}</h4>
</div>

<div class="modal-body">
    <div class="row">
        <div class="col-md-4">
            <a href="https://molgenis26.gcc.rug.nl/downloads/gavin/variantclassificationwebtool/plots_r0.3/${entity.get('Gene_Name')}.png" target="_blank"><img src="https://molgenis26.gcc.rug.nl/downloads/gavin/variantclassificationwebtool/plots_r0.3/${entity.get('Gene_Name')}.png" height="675" width="1200"/></a>
		</div>
    </div>
</div>

<#-- modal footer -->
<div class="modal-footer">
    <button type="button" class="btn btn-default" data-dismiss="modal">close</button>
</div>
