import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
    base: '/BiodiversityObservationNetworks.jl/previews/PR83/',// TODO: replace this in makedocs!
    title: 'BiodiversityObservationNetworks.jl',
    description: 'Documentation for BiodiversityObservationNetworks.jl',
    cleanUrls: true,
    outDir: '../1', // This is required for MarkdownVitepress to work correctly...
    head: [['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }]],
    ignoreDeadLinks: true,
    markdown: {
        math: true,
        config(md) {
            md.use(tabsMarkdownPlugin),
                md.use(mathjax3),
                md.use(footnote)
        },
        theme: {
            light: "min-light",
            dark: "min-dark"
        }
    },
    themeConfig: {
        outline: 'deep',
        siteTitle: 'BONs.jl',
        docFooter: {
            next: false,
            prev: false
        },
        search: {
            provider: 'local',
            options: {
                detailedView: true
            }
        },
        logo: { src: '/logo.png', width: 24, height: 24},
        nav: [
            {   
                text: "Tutorials", 
                items: [
                    { text: "Getting Started with BONs.jl", link: "/tutorials/gettingstarted.md"},
                    { text: "Building Multistage Samplers", link: "/tutorials/multistage.md"},
                    { text: "Measuring Spatial Balance", link: "/tutorials/spatialbalance.md"},
                    { text: "Measuring Environmental Representativeness", link: "/tutorials/envdistance.md"},
                    { text: "Designing a National BON", link: "/tutorials/canbon.md"},
                ],
            },
            { 
                text: "How-to", 
                link: "/howto",
            },
            { 
                text: "Reference", 
                items: [
                    { text: "Samplers", link: "/reference/samplers" },
                    { text: "Utilities", link: "/reference/utilities" },
                    { text: "API Reference", link: "/reference/api" },
                ]
            }
        ],
        sidebar: {
            "/tutorials/": [
                {
                    text: "Tutorials",
                    items: [
                        { text: "Getting Started with BONs.jl", link: "/tutorials/gettingstarted.md"},
                        { text: "Building Multistage Samplers", link: "/tutorials/multistage.md"},
                        { text: "Measuring Spatial Balance", link: "/tutorials/spatialbalance.md"},
                        { text: "Measuring Environmental Representativeness", link: "/tutorials/envdistance.md"},
                        { text: "Designing a National BON", link: "/tutorials/canbon.md"},
                    ],
                }
            ],
            "/howto/": [
                {
                    text: "How-to",
                    items: [],
                }
            ],
            "/reference/": [
                {
                    text: "Samplers",
                    collapsed: false,
                    items: [
                        { text: "Simple Random Sampling", link: "/reference/samplers/simplerandom" },
                        { text: "Balanced Acceptance Sampling", link: "/reference/samplers/balancedacceptance" },
                        { text: "Generalized Random Tessellated Stratified Sampling", link: "/reference/samplers/grts" },
                        { text: "Cube Sampling", link: "/reference/samplers/cube" },
                        { text: "Adaptive Hotspot Detection", link: "/reference/samplers/adaptivehotspot" },
                    ]
                },
                {
                    text: "Utilities",
                    collapsed: false,
                    items: [
                        { text: "Spatial Balance", link: "/reference/utilities/spatialbalance.md"},
                        { text: "Environmental Distance", link: "/reference/utilities/envdistance.md"},
                        { text: "Climate Rarity", link: "/reference/utilities/rarity.md"},
                        { text: "Climate Velocity", link: "/reference/utilities/velocity.md"},
                    ]
                },
                {
                    text: "Full API",
                    link: "/reference/api"
                }
            ]
        },
        editLink: { pattern: "https://https://github.com/PoisotLab/BiodiversityObservationNetworks.jl/edit/main/docs/src/:path" },
        socialLinks: [
            { icon: 'github', link: 'https://github.com/PoisotLab/BiodiversityObservationNetworks.jl' }
        ],
        footer: {
            message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> by the <a href="https://poisotlab.io/" target="_blank">Computational Ecology Research Group</a><br>',
            copyright: `This documentation is released under the <a href="https://creativecommons.org/licenses/by/4.0/deed.en" target="_blank">CC-BY 4.0 licence</a> - ${new Date().getUTCFullYear()}`
        }
    }
})
