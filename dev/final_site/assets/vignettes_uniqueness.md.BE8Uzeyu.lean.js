import{_ as a,c as s,a5 as n,o as t}from"./chunks/framework.BdiZ026L.js";const h=JSON.parse('{"title":"Selecting environmentally unique locations","description":"","frontmatter":{},"headers":[],"relativePath":"vignettes/uniqueness.md","filePath":"vignettes/uniqueness.md","lastUpdated":null}'),i={name:"vignettes/uniqueness.md"};function l(o,e,p,c,d,r){return t(),s("div",null,e[0]||(e[0]=[n(`<h1 id="Selecting-environmentally-unique-locations" tabindex="-1">Selecting environmentally unique locations <a class="header-anchor" href="#Selecting-environmentally-unique-locations" aria-label="Permalink to &quot;Selecting environmentally unique locations {#Selecting-environmentally-unique-locations}&quot;">​</a></h1><p>For some applications, we want to sample a set of locations that cover a broad range of values in environment space. Another way to rephrase this problem is to say we want to find the set of points with the <em>least</em> covariance in their environmental values.</p><p>To do this, we use a <code>BONRefiner</code> called <code>Uniqueness</code>. We&#39;ll start by loading the required packages.</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>using BiodiversityObservationNetworks</span></span>
<span class="line"><span>using SpeciesDistributionToolkit</span></span>
<span class="line"><span>using StatsBase</span></span>
<span class="line"><span>using NeutralLandscapes</span></span>
<span class="line"><span>using CairoMakie</span></span></code></pre></div><div class="warning custom-block"><p class="custom-block-title">Consider setting your SDMLAYERS_PATH</p><p>When accessing data using <code>SimpleSDMDatasets.jl</code>, it is best to set the <code>SDM_LAYERSPATH</code> environmental variable to tell <code>SimpleSDMDatasets.jl</code> where to download data. This can be done by setting <code>ENV[&quot;SDMLAYERS_PATH&quot;] = &quot;/home/user/Data/&quot;</code> or similar in the <code>~/.julia/etc/julia/startup.jl</code> file. (Note this will be different depending on where <code>julia</code> is installed.)</p></div><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>bbox = (left=-83.0, bottom=46.4, right=-55.2, top=63.7);</span></span>
<span class="line"><span>temp, precip, elevation = </span></span>
<span class="line"><span>    convert(Float32, SimpleSDMPredictor(RasterData(WorldClim2, AverageTemperature); bbox...)),</span></span>
<span class="line"><span>    convert(Float32, SimpleSDMPredictor(RasterData(WorldClim2, Precipitation); bbox...)),</span></span>
<span class="line"><span>    convert(Float32, SimpleSDMPredictor(RasterData(WorldClim2, Elevation); bbox...));</span></span></code></pre></div><p>Now we&#39;ll use the <code>stack</code> function to combine our four environmental layers into a single, 3-dimensional array, which we&#39;ll pass to our <code>Uniqueness</code> refiner.</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>layers = BiodiversityObservationNetworks.stack([temp,precip,elevation]);</span></span></code></pre></div><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>uncert = rand(MidpointDisplacement(0.8), size(temp), mask=temp);</span></span>
<span class="line"><span>heatmap(uncert)</span></span></code></pre></div><p>Now we&#39;ll get a set of candidate points from a BalancedAcceptance seeder that has no bias toward higher uncertainty values.</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>candpts = seed(BalancedAcceptance(numsites=100));</span></span></code></pre></div><p>Now we&#39;ll <code>refine</code> our <code>100</code> candidate points down to the 30 most environmentally unique.</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>finalpts = refine(candpts, Uniqueness(;numsites=30, layers=layers))</span></span>
<span class="line"><span>heatmap(uncert)</span></span>
<span class="line"><span>scatter!([p[1] for p in candpts], [p[2] for p in candpts], color=:white)</span></span>
<span class="line"><span>scatter!([p[1] for p in finalpts], [p[2] for p in finalpts], color=:dodgerblue, msc=:white)</span></span>
<span class="line"><span>current_figure()</span></span></code></pre></div>`,13)]))}const g=a(i,[["render",l]]);export{h as __pageData,g as default};