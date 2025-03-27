import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
    base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
    title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    cleanUrls: true,
    outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
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
        logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
        nav: [
            {   
                text: "Tutorials", 
                items: [
                    { text: "Getting Started with BONs.jl", link: "/tutorials/gettingstarted.md"},
                    { text: "Building Multistage Samplers", link: "/tutorials/multistage.md"},
                    { text: "Measuring Spatial Balance", link: "/tutorials/gettingstarted.md"},
                    { text: "Measuring Environmental Representativeness", link: "/tutorials/spatialbalance.md"},
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
                        { text: "Measuring Spatial Balance", link: "/tutorials/gettingstarted.md"},
                        { text: "Measuring Environmental Representativeness", link: "/tutorials/spatialbalance.md"},
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
        editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
        socialLinks: [
            { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
        ],
        footer: {
            message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> by the <a href="https://poisotlab.io/" target="_blank">Computational Ecology Research Group</a><br>',
            copyright: `This documentation is released under the <a href="https://creativecommons.org/licenses/by/4.0/deed.en" target="_blank">CC-BY 4.0 licence</a> - ${new Date().getUTCFullYear()}`
        }
    }
})
