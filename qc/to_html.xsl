<xsl:stylesheet id="to_html" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:ns="http://www.prime-xs.eu/ms/qcml" xmlns="">
    <xsl:template match="/">
        <html style="font-family:Arial;">
            <head>
                <link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css"/>
                <style>
                    table {
                        font-size: 10px;
                    }
                    th {
                        background-color: #CCCCCC; 
                    }
                    .table th {
                        text-align: center;   
                    }
                    .table td {
                        text-align: center;   
                    }
                    .border {
                        border-color: #000000;
                        border-width: 2px;
                        border-style: solid;
                        padding: 20px;
                        margin: 10px;
                    }
                    hr {
                        border-color: #000000;
                        margin: 25px 10px 0px 10px;
                    }
                </style>
            </head>
            <body>
                <div class="container text-center">
                    <div class="row">
                        <div class="col-md-12">
                            <h1>Analysis of mass spectrometry quality control metrics</h1>
                        </div>
                    </div>
                    <!-- overview -->
                    <div class="row border">
                        <div class="row">
                            <h2>Combined Analysis</h2>
                        </div>
                        <xsl:for-each select="ns:qcML/ns:setQuality">
                            <div class="row">
                                <div class="row">
                                    <h3>Visualization</h3>
                                </div>
                                <div class="row">
                                    <div class="col-md-4">
                                        <img class="img-responsive center-block"><xsl:attribute name="src">data:image/png;base64,
                                           <xsl:value-of select="ns:attachment[@ID='time']/ns:binary"/></xsl:attribute>
                                        </img>
                                        <b><xsl:value-of select="ns:attachment[@ID='time']/@name"/></b>
                                    </div>
                                    <div class="col-md-4">
                                       <img class="img-responsive center-block"><xsl:attribute name="src">data:image/png;base64,
                                           <xsl:value-of select="ns:attachment[@ID='PCA']/ns:binary"/></xsl:attribute>
                                       </img>
                                       <b><xsl:value-of select="ns:attachment[@ID='PCA']/@name"/></b>
                                    </div>
                                    <div class="col-md-4">
                                        <img class="img-responsive center-block"><xsl:attribute name="src">data:image/png;base64,
                                            <xsl:value-of select="ns:attachment[@ID='t-SNE']/ns:binary"/></xsl:attribute>
                                        </img>
                                        <b><xsl:value-of select="ns:attachment[@ID='t-SNE']/@name"/></b>
                                    </div>
                                </div>
                                <div class="row">
                                    <div class="col-md-3">
                                        <h3>Preprocessing</h3>
                                        <h4><xsl:value-of select="ns:attachment[@ID='var']/@name"/></h4>
                                        <p>
                                            <table class="table table-condensed">
                                                <tr>
                                                    <xsl:call-template name="table-header">
                                                        <xsl:with-param name="list"><xsl:value-of select="ns:attachment[@ID='var']/ns:table/ns:tableColumnTypes"/></xsl:with-param>
                                                    </xsl:call-template>
                                                </tr>
                                                <xsl:for-each select="ns:attachment[@ID='var']/ns:table/ns:tableRowValues">
                                                    <tr>
                                                        <xsl:call-template name="table-row">
                                                            <xsl:with-param name="list"><xsl:value-of select="." /></xsl:with-param>
                                                        </xsl:call-template>
                                                    </tr>
                                                </xsl:for-each>
                                            </table>
                                        </p>
                                        <p>
                                            <xsl:value-of select="ns:qualityParameter[@ID='VarianceThreshold']/@name"/> = <xsl:value-of select="ns:qualityParameter[@ID='VarianceThreshold']/@value"/>
                                        </p>
                                        <h4><xsl:value-of select="ns:attachment[@ID='corr']/@name"/></h4>
                                        <p>
                                            <table class="table table-condensed">
                                                <tr>
                                                    <xsl:call-template name="table-header">
                                                        <xsl:with-param name="list"><xsl:value-of select="ns:attachment[@ID='corr']/ns:table/ns:tableColumnTypes"/></xsl:with-param>
                                                    </xsl:call-template>
                                                </tr>
                                                <xsl:for-each select="ns:attachment[@ID='corr']/ns:table/ns:tableRowValues">
                                                    <tr>
                                                        <xsl:call-template name="table-row">
                                                            <xsl:with-param name="list"><xsl:value-of select="." /></xsl:with-param>
                                                        </xsl:call-template>
                                                    </tr>
                                                </xsl:for-each>
                                            </table>
                                        </p>
                                        <p>
                                            <xsl:value-of select="ns:qualityParameter[@ID='CorrelationThreshold']/@name"/> = <xsl:value-of select="format-number(ns:qualityParameter[@ID='CorrelationThreshold']/@value, '00.00%')"/>
                                        </p>
                                    </div>
                                    <div class="col-md-9">
                                        <div class="row">
                                            <h3>Outlier analysis</h3>
                                        </div>
                                        <div class="row">
                                            <div class="col-md-8">
                                                <figure>
                                                    <img class="img-responsive center-block"><xsl:attribute name="src">data:image/png;base64,
                                                        <xsl:value-of select="ns:attachment[@ID='OutlierScoreHistogram']/ns:binary"/></xsl:attribute>
                                                    </img>
                                                    <figcaption><b><xsl:value-of select="ns:qualityParameter[@ID='OutlierScoreThreshold']/@name"/> = <xsl:value-of select="format-number(ns:qualityParameter[@ID='OutlierScoreThreshold']/@value, '00.00%')"/><br/>
                                                        <xsl:value-of select="ns:qualityParameter[@ID='NrOutliers']/@name"/> = <xsl:value-of select="ns:qualityParameter[@ID='NrOutliers']/@value"/></b></figcaption>
                                                </figure>
                                            </div>
                                            <div class="col-md-4">
                                                <p>
                                                    <table class="table table-condensed">
                                                        <tr>
                                                            <xsl:call-template name="table-header">
                                                                <xsl:with-param name="list"><xsl:value-of select="ns:attachment[@ID='freq']/ns:table/ns:tableColumnTypes"/></xsl:with-param>
                                                            </xsl:call-template>
                                                        </tr>
                                                        <xsl:for-each select="ns:attachment[@ID='freq']/ns:table/ns:tableRowValues">
                                                            <tr>
                                                                <xsl:call-template name="table-row">
                                                                    <xsl:with-param name="list"><xsl:value-of select="." /></xsl:with-param>
                                                                </xsl:call-template>
                                                            </tr>
                                                        </xsl:for-each>
                                                    </table>
                                                </p>
                                                <p>
                                                    <xsl:value-of select="ns:qualityParameter[@ID='minsup']/@name"/> = <xsl:value-of select="ns:qualityParameter[@ID='minsup']/@value"/>
                                                </p>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </xsl:for-each>
                    </div>
                    <div class="row border">
                        <div class="row">
                            <h2>Individual Outliers</h2>
                        </div>
                        <xsl:for-each select="ns:qcML/ns:runQuality">
                            <div class="row">
                                <div class="row">
                                    <h3><xsl:value-of select="@ID"/></h3>
                                    <h5><xsl:value-of select="ns:qualityParameter[@name='Outlier score']/@name"/> = <xsl:value-of select="format-number(ns:qualityParameter[@name='Outlier score']/@value, '00.00%')"/></h5>
                                </div>
                                <div class="row">
                                    <div class="col-md-6">
                                        <img class="img-responsive center-block"><xsl:attribute name="src">data:image/png;base64,
                                            <xsl:value-of select="ns:attachment[@name='Feature importance']/ns:binary"/></xsl:attribute>
                                        </img>
                                        <b><xsl:value-of select="ns:attachment[@name='Feature importance']/@name"/></b>
                                    </div>
                                    <div class="col-md-6">
                                        <img class="img-responsive center-block"><xsl:attribute name="src">data:image/png;base64,
                                            <xsl:value-of select="ns:attachment[@name='Explanatory subspace']/ns:binary"/></xsl:attribute>
                                        </img>
                                        <b><xsl:value-of select="ns:attachment[@name='Explanatory subspace']/@name"/></b>                                    
                                    </div>
                                </div>
                            </div>
                            <xsl:choose>
                                <xsl:when test="position() != last()"><div class="row"><hr/></div></xsl:when>
                            </xsl:choose>
                        </xsl:for-each>
                    </div>
                </div>
            </body>
        </html>
    </xsl:template>
    
    <xsl:template name="table-header">
        <xsl:param name="list"/>
        <xsl:variable name="newlist" select="concat(normalize-space($list), ' ')"/>
        <xsl:variable name="first" select="substring-before($newlist, ' ')"/>
        <xsl:variable name="remaining" select="substring-after($newlist, ' ')"/>
        <th>
            <xsl:value-of select="$first"/>
        </th>
        <xsl:if test="$remaining">
            <xsl:call-template name="table-header">
                <xsl:with-param name="list" select="$remaining"/>
            </xsl:call-template>
        </xsl:if>
    </xsl:template>
    
    <xsl:template name="table-row">
        <xsl:param name="list"/>
        <xsl:variable name="newlist" select="concat(normalize-space($list), ' ')"/>
        <xsl:variable name="first" select="substring-before($newlist, ' ')"/>
        <xsl:variable name="remaining" select="substring-after($newlist, ' ')"/>
        <td>
            <xsl:value-of select="$first"/>
        </td>
        <xsl:if test="$remaining">
            <xsl:call-template name="table-row">
                <xsl:with-param name="list" select="$remaining"/>
            </xsl:call-template>
        </xsl:if>
    </xsl:template>
</xsl:stylesheet>
